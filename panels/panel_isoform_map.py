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
from matplotlib.patches import Rectangle, Patch
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from constants import FOLDER_PATH, BAND_RANGE, SAVE_PLOTS, BAND_FOLDER, FASTA_PATH, PEPTIDE_CSV_DIR
from constants import GENETIC_CODE, GEL_BANDS, AA_WEIGHTS, WATER_MASS
from utils import *
from ucsc_api import *
from PIL import Image

def analyze_isoforms(fasta_path, peptide_csv_dir, gene_name):
    fasta_seqs = parse_fasta(fasta_path, gene_name)
    isoform_to_peptides = find_peptide_to_isoform_map(peptide_csv_dir, gene_name)

    results = []

    for isoform_id in isoform_to_peptides:
        matching_keys = [k for k in fasta_seqs if isoform_id in k]
        if not matching_keys:
            continue
        key = matching_keys[0]
        seq = fasta_seqs[key]
        mw = calculate_mw(seq)
        bands = assign_band(mw)
        peptides = sorted(set(isoform_to_peptides[isoform_id]))
        results.append({
            "isoform_id": key,
            "mw_kda": round(mw, 2),
            "bands": bands,
            "sequence": seq,
            "peptides": peptides
        })

    results.sort(key=lambda r: r["mw_kda"])
    return results


def plot_isoform_maps(results, gene_name, output_dir):
    # Fetch canonical transcript and exon mapping info
    transcripts = fetch_all_transcripts(gene_name)
    tx = select_longest_cds_transcript(transcripts)
    if tx is None:
        raise ValueError(f"No transcript found for gene {gene_name}")
    canonical_protein, aa_ranges, _ = map_exons_to_aa(tx)

    # Strip terminal stop codon from canonical protein
    canonical_protein = canonical_protein.rstrip('*')

    # Sort aa_ranges to ensure N-terminal → C-terminal exon order
    aa_ranges = sorted(aa_ranges, key=lambda x: x[0])

    # Strip terminal stop codon from all isoform sequences
    for r in results:
        r["sequence"] = r["sequence"].rstrip('*')

    # Prepare colors
    colors = cm.get_cmap('tab20').colors

    def is_red(color, threshold=0.8):
        r, g, b = color[:3]
        return (r > threshold) and (g < 0.3) and (b < 0.3)

    filtered_colors = [c for c in colors if not is_red(c)]
    exon_color_map = {i: filtered_colors[i % len(filtered_colors)] for i in range(len(aa_ranges))}

    fig_height = 0.6 * len(results)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    y_offset = 0
    max_len = max(len(r["sequence"]) for r in results)
    used_exons = set()

    for r in results:
        iso_seq = r["sequence"]
        seq_len = len(iso_seq)
        iso_y = y_offset + 1

        # Draw dotted line using Line2D with short dash spacing
        from matplotlib.lines import Line2D
        line = Line2D(
            xdata=[0, seq_len],
            ydata=[iso_y, iso_y],
            color='black',
            linewidth=4
        )
        line.set_dashes([1, 0.5])  # dash, gap spaces
        ax.add_line(line)

        # Map exons from canonical protein to isoform
        for exon_idx, (start, end) in enumerate(aa_ranges):
            exon_seq = canonical_protein[start:end]
            pos_in_iso = iso_seq.find(exon_seq)
            if pos_in_iso != -1:
                ax.hlines(iso_y, pos_in_iso, pos_in_iso + len(exon_seq),
                          colors=exon_color_map[exon_idx], linewidth=4)
                used_exons.add(exon_idx)

        # Overlay peptide boxes colored by exon mapping
        for pep in r["peptides"]:
            pep_start = iso_seq.find(pep)
            if pep_start == -1:
                continue

            match = re.search(re.escape(pep), canonical_protein)
            if not match:
                ax.add_patch(Rectangle(
                    (pep_start, iso_y - 0.2), len(pep), 0.4,
                    facecolor='#FF0000', edgecolor='black', lw=0.7, hatch='////', zorder=3))
                continue

            start, end = match.start(), match.end()
            exon_indices = []
            for pos in range(start, end):
                for exon_idx, (exon_start, exon_end) in enumerate(aa_ranges):
                    if exon_start <= pos < exon_end:
                        exon_indices.append(exon_idx)
                        break
                else:
                    exon_indices.append(None)

            label_segments = []
            current_exon = exon_indices[0]
            current_segment = pep[0]
            for aa_idx, exon in zip(range(1, len(pep)), exon_indices[1:]):
                if exon == current_exon:
                    current_segment += pep[aa_idx]
                else:
                    label_segments.append((current_segment, current_exon))
                    current_segment = pep[aa_idx]
                    current_exon = exon
            label_segments.append((current_segment, current_exon))

            offset_pep = pep_start
            for seg, exon in label_segments:
                color = exon_color_map[exon] if exon is not None else 'black'
                if exon is not None:
                    used_exons.add(exon)
                ax.add_patch(Rectangle((offset_pep, iso_y - 0.2), len(seg), 0.4,
                                       facecolor=color, edgecolor='black', lw=0.5, zorder=3))
                offset_pep += len(seg)

        label = f"{r['isoform_id']}\n{r['mw_kda']} kDa, Bands: {', '.join(map(str, r['bands']))}"
        ax.text(-0.01 * max_len, iso_y, label, va="center", ha="right", fontsize=8)
        y_offset += 1.5

    ax.set_xlim(0, max_len)
    ax.set_ylim(0, y_offset + 1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"Detected Peptides in {gene_name} Isoforms (Exon-labeled)")

    # Legend: ensure N-terminal → C-terminal exon order
    exon_legend_patches = [
        Patch(facecolor=exon_color_map[i], label=f"Exon {i+1}")
        for i in range(len(aa_ranges)) if i in used_exons
    ]
    ax.legend(handles=exon_legend_patches, title="Exons", bbox_to_anchor=(1.02, 1),
              loc='upper left', borderaxespad=0.)

    plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{gene_name}_protein_isoform_map.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_file
