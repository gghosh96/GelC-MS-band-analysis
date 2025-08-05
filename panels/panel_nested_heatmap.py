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
from constants import FOLDER_PATH, BAND_RANGE, SAVE_PLOTS, BAND_FOLDER, FASTA_PATH, PEPTIDE_CSV_DIR
from constants import GENETIC_CODE, GEL_BANDS, AA_WEIGHTS, WATER_MASS
from utils import *
from ucsc_api import *
from PIL import Image
        
def get_peptide_exon_mapping(peptides, protein_seq, aa_to_exon):
    # Maps each peptide to the set of exons it spans based on amino acid positions
    mapping = {}
    for pep in peptides:
        match = re.search(re.escape(pep), protein_seq)
        if not match:
            continue
        start, end = match.start(), match.end()
        exons = sorted(set(aa_to_exon.get(i) for i in range(start, end) if i in aa_to_exon))
        mapping[pep] = exons
    return mapping

def plot_panel_b(protein_name, output_dir):
    # Load all band-specific peptide quantification files and sort them by band number
    band_files = sorted(
        [(int(re.search(r'Band_(\d+)', f).group(1)), f) for f in os.listdir(BAND_FOLDER) if f.endswith('.csv')],
        key=lambda x: x[0]
    )
    bands = [str(band_num) for band_num, _ in band_files]
    peptide_data = {}

    # Process each gel band to extract peptide intensity data for the specified protein
    for band_num, filename in band_files:
        path = os.path.join(BAND_FOLDER, filename)
        df = pd.read_csv(path)
        df['Protein'] = df['Protein'].astype(str)

        # Keep only rows matching the target protein (or isoforms)
        df = df[df['Protein'].str.split(';').apply(lambda prots: any(matches_gene_id(p.strip(), protein_name) for p in prots))]

        if df.empty:
            continue

        # Normalize intensities by max value across tumor/normal samples per peptide
        intensity_cols = [col for col in df.columns if '_T' in col or '_N' in col]
        max_intensity = df[intensity_cols].max(axis=1).replace(0, np.nan)

        # For each peptide, collect normalized tumor and normal values, and compute fold change
        for idx, row in df.iterrows():
            peptide = row['sequences']
            norm = row[intensity_cols] / max_intensity.loc[idx]
            T_vals = norm[[c for c in intensity_cols if '_T' in c]].dropna().tolist()
            N_vals = norm[[c for c in intensity_cols if '_N' in c]].dropna().tolist()
            if not T_vals or not N_vals:
                continue
            fc = np.log2(np.mean(T_vals) / np.mean(N_vals)) if np.mean(N_vals) > 0 else 0
            peptide_data.setdefault(peptide, {})[str(band_num)] = {'T': T_vals, 'N': N_vals, 'fold_change': fc}

    if not peptide_data:
        raise ValueError(f"No peptides found for {protein_name}")

    # Fetch transcript and exon structure from UCSC and map amino acids to exons
    tx = fetch_transcript_info(protein_name)
    protein_seq, aa_ranges, aa_to_exon = map_exons_to_aa(tx)

    # Order peptides based on their position in the protein sequence
    peptides = sorted(peptide_data, key=lambda p: protein_seq.find(p) if protein_seq.find(p) != -1 else float('inf'))

    # Map each peptide to the exon(s) it overlaps
    pep_to_exons = get_peptide_exon_mapping(peptides, protein_seq, aa_to_exon)
    exon_colors = plt.cm.tab20.colors

    # Prepare layout for nested heatmap grid: one row per peptide, one column per band
    n_rows = len(peptides)
    n_cols = len(bands)
    fig_width = max(12, n_cols * 1.3)
    fig_height = max(4, n_rows * 0.7)
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(n_rows, n_cols, figure=fig, wspace=0.05, hspace=0.05)

    # Get range of fold change values to normalize color scale across all bands/peptides
    all_fc = [data[band]['fold_change'] for data in peptide_data.values() for band in data]
    max_abs = max(abs(x) for x in all_fc) if all_fc else 1

    # For each peptide-band pair, create a mini heatmap cell with bars for tumor and normal
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
            ax.set_facecolor('white')
            ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=color, zorder=0))
            T_vals = peptide_data[pep][band]['T']
            N_vals = peptide_data[pep][band]['N']
            bar_width = 0.3
            x_N = list(range(len(N_vals)))
            x_T = [x + bar_width for x in x_N]
            ax.bar(x_N, N_vals, width=bar_width, color='black', edgecolor='none')
            ax.bar(x_T, T_vals, width=bar_width, color='lightcoral', edgecolor='none')
            combined = T_vals + N_vals
            if combined:
                ax.set_ylim(min(combined) * 0.9, max(combined) * 1.1)
            ax.set_aspect('auto')

    # Annotate each row with the peptide sequence, using exon-colored segmentation if available
    label_x = -0.12
    char_width = 0.0065
    gap = 0.0005
    for i, pep in enumerate(peptides):
        y_pos = (n_rows - i - 0.5) / n_rows
        match = re.search(re.escape(pep), protein_seq)
        if not match:
            fig.text(label_x, y_pos, pep, va='center', ha='left',
                     fontsize=12, fontweight='bold', fontfamily='monospace', color='black')
            continue
        start, end = match.start(), match.end()
        exon_idxs = [aa_to_exon.get(pos, None) for pos in range(start, end)]

        # Split peptide into segments based on exon boundaries to color accordingly
        segs, current_exon, current_seg = [], exon_idxs[0], pep[0]
        for aa_idx, exon in zip(range(1, len(pep)), exon_idxs[1:]):
            if exon == current_exon:
                current_seg += pep[aa_idx]
            else:
                segs.append((current_seg, current_exon))
                current_seg, current_exon = pep[aa_idx], exon
        segs.append((current_seg, current_exon))

        # Render peptide segments with exon-specific colors
        offset_x = label_x
        for seg, exon in segs:
            color = exon_colors[exon % len(exon_colors)] if exon is not None else 'gray'
            fig.text(offset_x, y_pos, seg, va='center', ha='left',
                     fontsize=16, fontweight='bold', fontfamily='monospace', color=color)
            offset_x += len(seg) * char_width + gap

    # Label each gel band column
    for j, band in enumerate(bands):
        fig.text((j + 0.5) / n_cols, -0.01, f'Band {band}', va='top', ha='center', fontsize=12)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # Create colorbar legend to represent log2 fold-change values
    norm = mcolors.Normalize(vmin=-max_abs, vmax=max_abs)
    sm = cm.ScalarMappable(cmap='BrBG', norm=norm)
    sm.set_array([])

    # Align colorbar height with peptide heatmap
    fig_width, fig_height = fig.get_size_inches()
    row_height_inch = fig_height / n_rows
    heatmap_height_inch = n_rows * row_height_inch
    cbar_height_frac = heatmap_height_inch / fig_height
    cbar_bottom = (1 - cbar_height_frac) / 2
    cbar_ax = fig.add_axes([1.01, cbar_bottom, 0.015, cbar_height_frac])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Log2 Fold Change (Tumor / Normal)', fontsize=12)

    # Add main figure title with protein name
    title_offset = (fig_height + 0.4) / fig_height
    fig.suptitle(f"{protein_name} Peptide Heatmap (Tumor vs Normal)",
                 fontsize=18, y=title_offset)

    # Save the final heatmap figure
    output_file = output_dir / f"{protein_name}_peptide_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_file
