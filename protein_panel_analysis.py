import os
import re
import sys
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
from panels.panel_boxplot import *
from panels.panel_isoform_map import *
from panels.panel_nested_heatmap import *
from panels.panel_protein_map import *
from PIL import Image
from panels.panel_clustered_isoform_map import (
    analyze_isoforms_clustered,
    plot_clustered_isoform_maps_with_sas
)
from panels.panel_log_boxplot import plot_log_protein_boxplot   # <-- NEW IMPORT


def analyze_all_panels(gene_name):
    output_dir = Path("data") / gene_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------
    # protein intensity boxplot
    # -----------------------------
    print(f"Loading data and plotting Panel A for {gene_name}...")
    df = load_diffacto_data(gene_name)
    panelA_file = plot_panel_a_boxplot(df, gene_name, output_dir)

    # -----------------------------
    # NEW: log intensity boxplot + t-test results
    # -----------------------------
    print(f"Plotting log-intensity boxplot and calculating t-tests for {gene_name}...")
    try:
        log_boxplot_file, log_ttest_file = plot_log_protein_boxplot(df, gene_name, output_dir)
    except ValueError as e:
        print(f"Skipping log-intensity plot: {e}")
        log_boxplot_file, log_ttest_file = None, None

    # -----------------------------
    # protein exon mapping plot
    # -----------------------------
    print(f"Fetching transcript and plotting exon-colored protein for {gene_name}...")
    transcripts = fetch_all_transcripts(gene_name)
    tx = select_longest_cds_transcript(transcripts)
    protein_seq, aa_ranges, aa_to_exon = map_exons_to_aa(tx)
    panelProtein_file = plot_protein(protein_seq, aa_ranges, gene_name, output_dir=output_dir)

    # -----------------------------
    # nested peptide heatmap
    # -----------------------------  
    print(f"Plotting Panel B heatmap for {gene_name}...")
    panelB_file = plot_panel_b(gene_name, output_dir)

    # -----------------------------
    # Original isoform map
    # -----------------------------
    print(f"Analyzing isoforms and plotting exon-colored maps for {gene_name}...")
    isoform_results = analyze_isoforms(FASTA_PATH, PEPTIDE_CSV_DIR, gene_name)
    isoform_map_file = plot_isoform_maps(isoform_results, gene_name, output_dir)

    # -----------------------------
    # Clustered isoform map (with truncations/substitutions)
    # -----------------------------
    print(f"Analyzing isoforms (clustered) and plotting clustered peptide maps for {gene_name}...")
    clustered_results = analyze_isoforms_clustered(FASTA_PATH, PEPTIDE_CSV_DIR, gene_name)
    clustered_file, cluster_txt_file = plot_clustered_isoform_maps_with_sas(
        clustered_results, gene_name, output_dir
    )

    print("All panels generated and saved individually.")

    return {
        'panelA': panelA_file,
        'panelLogBoxplot': log_boxplot_file,      # <-- NEW OUTPUT
        'log_ttest_csv': log_ttest_file,          # <-- NEW OUTPUT
        'panelProtein': panelProtein_file,
        'panelB': panelB_file,
        'isoform_map': isoform_map_file,
        'clustered_isoform_map': clustered_file,
        'cluster_txt': cluster_txt_file
    }


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python main.py <GENE_NAME>")
        sys.exit(1)

    gene = sys.argv[1]
    panel_files = analyze_all_panels(gene)

    for panel_name, filepath in panel_files.items():
        print(f"{panel_name} saved at: {filepath}")
        if filepath is not None and panel_name not in ['cluster_txt', 'log_ttest_csv']:
            img = Image.open(filepath)
            img.show()