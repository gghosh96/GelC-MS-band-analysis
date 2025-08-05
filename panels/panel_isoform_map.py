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

def analyze_isoforms(fasta_path, peptide_csv_dir, gene_name):
    # Parse all protein sequences for a given gene from the FASTA file
    fasta_seqs = parse_fasta(fasta_path, gene_name)

    # Create a mapping of isoform IDs to peptides detected from peptide search results
    isoform_to_peptides = find_peptide_to_isoform_map(peptide_csv_dir, gene_name)

    results = []

    for isoform_id in isoform_to_peptides:
        # Match isoform ID with corresponding FASTA sequence entry
        matching_keys = [k for k in fasta_seqs if isoform_id in k]
        for key in matching_keys:
            seq = fasta_seqs[key]

            # Compute molecular weight and assign band(s) based on gel separation predictions
            mw = calculate_mw(seq)
            bands = assign_band(mw)

            # Collect peptide list for this isoform
            peptides = sorted(isoform_to_peptides[isoform_id])

            results.append({
                "isoform_id": key,
                "mw_kda": round(mw, 2),
                "bands": bands,
                "sequence": seq,
                "peptides": peptides
            })

    # Sort isoforms by molecular weight for visual consistency
    results.sort(key=lambda r: r["mw_kda"])
    return results


def plot_isoform_maps(results, gene_name, output_dir):
    # Fetch canonical transcript and exon mapping information from UCSC or equivalent source
    tx = fetch_transcript_info(gene_name)
    canonical_protein, _, aa_to_exon = map_exons_to_aa(tx)

    exon_colors = plt.cm.tab20.colors  # color map for exon-based highlighting

    fig_height = 0.6 * len(results)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    y_offset = 0
    max_len = max(len(r["sequence"]) for r in results)
    used_exons = set()

    for r in results:
        # Draw base line for isoform protein length
        seq_len = len(r["sequence"])
        iso_y = y_offset + 1
        ax.hlines(iso_y, 0, seq_len, color="black", linewidth=2)

        for pep in r["peptides"]:
            # Locate peptide in the isoform sequence
            pep_start = r["sequence"].find(pep)
            if pep_start == -1:
                continue

            # Map peptide to canonical protein to determine exon origin
            match = re.search(re.escape(pep), canonical_protein)
            if not match:
                # Peptides not mapping to canonical sequence are shown in gray
                ax.add_patch(plt.Rectangle((pep_start, iso_y - 0.2), len(pep), 0.4,
                                           facecolor='gray', edgecolor='black', lw=0.5))
                continue

            # Get exon indices for each residue in the peptide
            start, end = match.start(), match.end()
            exon_indices = [aa_to_exon.get(pos, None) for pos in range(start, end)]

            # Split peptide into exon-specific segments for coloring
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

            # Draw colored blocks for each exon segment in the peptide
            offset = pep_start
            for seg, exon in label_segments:
                color = exon_colors[exon % len(exon_colors)] if exon is not None else 'gray'
                if exon is not None:
                    used_exons.add(exon)
                ax.add_patch(plt.Rectangle((offset, iso_y - 0.2), len(seg), 0.4,
                                           facecolor=color, edgecolor='black', lw=0.5))
                offset += len(seg)

        # Annotate isoform ID, molecular weight, and predicted gel bands
        label = f"{r['isoform_id']}\n{r['mw_kda']} kDa, Bands: {', '.join(map(str, r['bands']))}"
        ax.text(-0.01 * max_len, iso_y, label, va="center", ha="right", fontsize=8)
        y_offset += 1.5

    ax.set_xlim(0, max_len)
    ax.set_ylim(0, y_offset + 1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"Detected Peptides in {gene_name} Isoforms (Exon-labeled)")

    # Create a legend showing colors assigned to each exon used
    from matplotlib.patches import Patch
    exon_legend_patches = [
        Patch(facecolor=exon_colors[exon % len(exon_colors)], label=f"Exon {exon}")
        for exon in sorted(used_exons)
    ]
    ax.legend(handles=exon_legend_patches, title="Exons", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

    plt.tight_layout()
    output_file = output_dir / f"{gene_name}_protein_isoform_map.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_file
