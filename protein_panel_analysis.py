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

def analyze_all_panels(gene_name):
    output_dir = Path("data") / gene_name
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data and plotting Panel A for {gene_name}...")
    df = load_diffacto_data(gene_name)
    panelA_file = plot_panel_a_boxplot(df, gene_name, output_dir)

    print(f"Fetching transcript and plotting exon-colored protein for {gene_name}...")
    tx = fetch_transcript_info(gene_name)
    protein_seq, aa_ranges, aa_to_exon = map_exons_to_aa(tx)
    panelProtein_file = plot_protein(protein_seq, aa_ranges, gene_name, output_dir=output_dir)

    print(f"Plotting Panel B heatmap for {gene_name}...")
    panelB_file = plot_panel_b(gene_name, output_dir)

    print(f"Analyzing isoforms and plotting isoform peptide maps for {gene_name}...")
    isoform_results = analyze_isoforms(FASTA_PATH, PEPTIDE_CSV_DIR, gene_name)
    isoform_file = plot_isoform_maps(isoform_results, gene_name, output_dir)

    print("All panels generated and saved individually.")
    
    return {
        'panelA': panelA_file,
        'panelProtein': panelProtein_file,
        'panelB': panelB_file,
        'isoform': isoform_file
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main.py <GENE_NAME>")
        sys.exit(1)

    gene = sys.argv[1]
    panel_files = analyze_all_panels(gene)

    for panel_name, filepath in panel_files.items():
        print(f"Displaying {panel_name} saved at: {filepath}")
        img = Image.open(filepath)
        img.show()