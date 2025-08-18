import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import mysql.connector
from pathlib import Path
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.patches import FancyArrowPatch
from constants import FOLDER_PATH, BAND_RANGE, SAVE_PLOTS, BAND_FOLDER, FASTA_PATH, PEPTIDE_CSV_DIR
from constants import GENETIC_CODE, GEL_BANDS, AA_WEIGHTS, WATER_MASS
from utils import *
from ucsc_api import *
from PIL import Image

def get_filtered_colors():
    colors = cm.get_cmap('tab20').colors

    def is_red(color, threshold=0.8):
        r, g, b = color[:3]
        return (r > threshold) and (g < 0.3) and (b < 0.3)

    filtered_colors = [c for c in colors if not is_red(c)]
    return filtered_colors

def plot_panel_b(protein_name, output_dir):
    # Load all band-specific peptide quantification files and sort them by band number
    band_files = sorted(
        [(int(re.search(r'Band_(\d+)', f).group(1)), f) for f in os.listdir(BAND_FOLDER) if f.endswith('.csv')],
        key=lambda x: x[0]
    )
    bands = [str(band_num) for band_num, _ in band_files]
    peptide_data = {}

    # First pass: load all peptide data for the protein
    all_data = []
    for band_num, filename in band_files:
        path = os.path.join(BAND_FOLDER, filename)
        df = pd.read_csv(path)
        df['Protein'] = df['Protein'].astype(str)
        df = df[df['Protein'].str.split(';').apply(lambda prots: any(matches_gene_id(p.strip(), protein_name) for p in prots))]
        if df.empty:
            continue
        df['Band'] = str(band_num)
        all_data.append(df)
    if not all_data:
        raise ValueError(f"No peptides found for {protein_name}")
    all_data = pd.concat(all_data, ignore_index=True)

    # Identify all intensity columns (tumor and normal)
    intensity_cols = [col for col in all_data.columns if '_T' in col or '_N' in col]

    # Normalize original intensities per patient across all peptides (max-normalization) - for fold change calculation
    patient_ids = sorted(set(col.split('_')[0] for col in intensity_cols))
    for pid in patient_ids:
        patient_cols = [c for c in intensity_cols if c.startswith(pid + '_')]
        max_val = all_data[patient_cols].to_numpy().max()
        if max_val > 0:
            all_data[patient_cols] = all_data[patient_cols] / max_val

    # Create fraction-normalized copy for bar graph proportions (sum to 1 per peptide per patient)
    fraction_data = all_data.copy()
    for pid in patient_ids:
        patient_cols = [c for c in intensity_cols if c.startswith(pid + '_')]
        for peptide, group_df in fraction_data.groupby('sequences'):
            vals = group_df[patient_cols].to_numpy()
            total_sum = vals.sum()
            if total_sum > 0:
                fraction_data.loc[group_df.index, patient_cols] = group_df[patient_cols] / total_sum

    # Build peptide_data combining fold change from all_data and fraction bars from fraction_data
    for _, row in all_data.iterrows():
        peptide = row['sequences']
        band_num = row['Band']

        # Fold change calculation with max-normalized intensities (all_data)
        T_vals = [row[c] for c in intensity_cols if c.startswith(tuple(patient_ids)) and '_T' in c]
        N_vals = [row[c] for c in intensity_cols if c.startswith(tuple(patient_ids)) and '_N' in c]
        T_vals = [v for v in T_vals if pd.notna(v)]
        N_vals = [v for v in N_vals if pd.notna(v)]
        if not T_vals or not N_vals:
            continue
        fc = np.log2(np.mean(T_vals) / np.mean(N_vals)) if np.mean(N_vals) > 0 else 0

        # Fraction values for bar graphs (fraction_data)
        T_frac = [fraction_data.loc[row.name, c] for c in intensity_cols if c.startswith(tuple(patient_ids)) and '_T' in c]
        N_frac = [fraction_data.loc[row.name, c] for c in intensity_cols if c.startswith(tuple(patient_ids)) and '_N' in c]

        peptide_data.setdefault(peptide, {})[band_num] = {'T': T_frac, 'N': N_frac, 'fold_change': fc}

    # Fetch transcript and exon mapping
    transcripts = fetch_all_transcripts(protein_name)
    tx = select_longest_cds_transcript(transcripts)
    if tx is None:
        raise ValueError(f"No transcript found for gene {protein_name}")
    protein_seq, aa_ranges, aa_to_exon = map_exons_to_aa(tx)

    # Order peptides by position in protein
    peptides = sorted(peptide_data, key=lambda p: protein_seq.find(p) if protein_seq.find(p) != -1 else float('inf'))

    # --- FIXED: Setup exon colors consistent with isoform map and strand orientation ---
    filtered_colors = get_filtered_colors()
    n_exons = len(aa_ranges)

    # Get strand safely from tx (handle dict or object)
    strand = getattr(tx, "strand", None)
    if strand is None and isinstance(tx, dict):
        strand = tx.get("strand", "+")
    if strand not in ["+", "-"]:
        strand = "+"

    exon_color_map = {}
    for exon_idx in range(n_exons):
        palette_idx = exon_idx if strand == '+' else (n_exons - 1 - exon_idx)
        exon_color_map[exon_idx] = filtered_colors[palette_idx % len(filtered_colors)]
    # -------------------------------------------------------------------

    def get_exon_index_for_position(pos, aa_ranges):
        # pos is 0-based (Python string index), aa_ranges are 1-based inclusive
        pos_1based = pos - 1
        for exon_idx, (exon_start, exon_end) in enumerate(aa_ranges):
            if exon_start <= pos_1based <= exon_end:  # inclusive on both ends
                return exon_idx
        return None

    # Figure setup
    n_rows, n_cols = len(peptides), len(bands)
    fig_width, fig_height = max(12, n_cols * 1.3), max(4, n_rows * 0.7)
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(n_rows, n_cols, figure=fig, wspace=0.05, hspace=0.05)

    # Calculate max absolute fold change for color scaling
    all_fc = [data[band]['fold_change'] for data in peptide_data.values() for band in data]
    max_abs = max(abs(x) for x in all_fc) if all_fc else 1

    # Calculate max total stacked bar height per band across all peptides for scaling
    max_total_per_band = {}
    for band in bands:
        totals = []
        for pep in peptides:
            if band in peptide_data.get(pep, {}):
                T_vals = peptide_data[pep][band]['T']
                N_vals = peptide_data[pep][band]['N']
                total = sum(t + n for t, n in zip(T_vals, N_vals))
                totals.append(total)
        max_total_per_band[band] = max(totals) if totals else 1.0  # avoid zero division

    # Plot mini heatmaps with stacked bars scaled per band max
    for i, pep in enumerate(peptides):
        for j, band in enumerate(bands):
            ax = fig.add_subplot(gs[i, j])
            ax.set_xticks([])
            ax.set_yticks([])

            if band not in peptide_data[pep]:
                ax.set_facecolor('white')
                continue

            fc = peptide_data[pep][band]['fold_change']
            color = plt.cm.BrBG((fc + max_abs) / (2 * max_abs))
            ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=color, zorder=0))

            T_vals = peptide_data[pep][band]['T']
            N_vals = peptide_data[pep][band]['N']

            total_vals = [t + n for t, n in zip(T_vals, N_vals)]

            # Scale stacked bars to max total per band
            max_total = max_total_per_band.get(band, 1.0)
            if max_total > 0:
                T_vals = [v / max_total for v in T_vals]
                N_vals = [v / max_total for v in N_vals]
                total_vals = [v / max_total for v in total_vals]

            bar_width = 0.6
            x = list(range(len(total_vals)))

            ax.bar(x, N_vals, width=bar_width, color='black', edgecolor='none', label='Normal' if (i == 0 and j == 0) else "")
            ax.bar(x, T_vals, width=bar_width, bottom=N_vals, color='lightcoral', edgecolor='none', label='Tumor' if (i == 0 and j == 0) else "")

            if total_vals:
                ax.set_ylim(0, max(total_vals) * 1.1)
            ax.set_aspect('auto')

    # Add legend once
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')

    char_width, gap = 0.0065, 0.0002

    # Calculate longest peptide total width
    longest_pep = max(peptides, key=len)
    longest_len = len(longest_pep)
    longest_total_width = longest_len * char_width + (longest_len - 1) * gap

    # Desired right edge position of peptide labels (0.5 left of vertical line at -0.005)
    label_right_edge = -0.005 - 0.005

    # Peptide label coloring with consistent exon colors
    for i, pep in enumerate(peptides):
        y_pos = (n_rows - i - 0.5) / n_rows
        match = re.search(re.escape(pep), protein_seq)
        if not match:
            pep_len = len(pep)
            total_width = pep_len * char_width + (pep_len - 1) * gap
            start_x = label_right_edge - total_width

            # Draw peptide in stop-sign red
            fig.text(start_x, y_pos, pep, va='center', ha='left',
                    fontsize=16, fontweight='bold', fontfamily='monospace',
                    color='#FF0000')

            # Underline (drawn as a separate line)
            underline_y = y_pos - 0.012  # tweak vertical position if needed
            underline_x_start = start_x
            underline_x_end = start_x + pep_len * char_width + (pep_len - 1) * gap
            fig.add_artist(plt.Line2D(
                [underline_x_start, underline_x_end],
                [underline_y, underline_y],
                transform=fig.transFigure,
                color='#FF0000', linewidth=1.5
            ))

            # Add arrow pointing to peptide
            arrow = FancyArrowPatch(
                posA=(start_x - 0.03, y_pos),   # arrow tail
                posB=(start_x - 0.004, y_pos),   # arrow head
                transform=fig.transFigure,      # figure fraction coords
                arrowstyle='->',
                mutation_scale=15,
                color='#FF0000',
                linewidth=1.5
            )
            fig.patches.append(arrow)
            continue

        start, end = match.start(), match.end()

        exon_idxs = [get_exon_index_for_position(pos + 1, aa_ranges) 
             for pos in range(start, end)]

        segs, current_exon, current_seg = [], exon_idxs[0], pep[0]
        for aa_idx, exon in zip(range(1, len(pep)), exon_idxs[1:]):
            if exon == current_exon:
                current_seg += pep[aa_idx]
            else:
                segs.append((current_seg, current_exon))
                current_seg, current_exon = pep[aa_idx], exon
        segs.append((current_seg, current_exon))

        pep_len = len(pep)
        total_width = pep_len * char_width + (pep_len - 1) * gap
        start_x = label_right_edge - total_width

        offset_x = start_x
        for seg, exon in segs:
            color = '#FF0000' if exon is None else exon_color_map.get(exon, 'black')  # dark red if no exon
            fig.text(offset_x, y_pos, seg, va='center', ha='left',
                     fontsize=16, fontweight='bold', fontfamily='monospace', color=color)
            offset_x += len(seg) * char_width + gap

    # Add vertical boundary line between labels and heatmap
    fig.add_artist(plt.Line2D(
        [-0.005, -0.005],
        [0, 1],
        transform=fig.transFigure,
        color='black', linewidth=1
    ))

    # Band labels
    for j, band in enumerate(bands):
        fig.text((j + 0.5) / n_cols, -0.01, f'Band {band}', va='top', ha='center', fontsize=12)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # Colorbar
    norm = mcolors.Normalize(vmin=-max_abs, vmax=max_abs)
    sm = cm.ScalarMappable(cmap='BrBG', norm=norm)
    sm.set_array([])
    fig_width, fig_height = fig.get_size_inches()
    row_height_inch = fig_height / n_rows
    heatmap_height_inch = n_rows * row_height_inch
    cbar_height_frac = heatmap_height_inch / fig_height
    cbar_bottom = (1 - cbar_height_frac) / 2
    cbar_ax = fig.add_axes([1.01, cbar_bottom, 0.015, cbar_height_frac])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Log2 Fold Change (Tumor / Normal)', fontsize=12)

    # Title
    title_offset = (fig_height + 0.4) / fig_height
    fig.suptitle(f"{protein_name} Peptide Heatmap (Tumor vs Normal)", fontsize=18, y=title_offset)

    # Save
    output_file = output_dir / f"{protein_name}_peptide_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_file
