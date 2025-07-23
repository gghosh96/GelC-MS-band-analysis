import os
import re
import numpy as np
import pandas as pd
import requests
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from constants import FOLDER_PATH, BAND_RANGE
from constants import GENETIC_CODE, GEL_BANDS, AA_WEIGHTS, WATER_MASS

def matches_gene_id(fasta_id, gene_name):
    pattern = re.compile(rf'^{re.escape(gene_name)}\.', re.IGNORECASE)
    return pattern.match(fasta_id) is not None


def load_diffacto_data(protein_name):
    all_data = []
    for band in BAND_RANGE:
        file_name = f"Band_{band}_proteins.csv"
        file_path = os.path.join(FOLDER_PATH, file_name)
        if not os.path.exists(file_path):
            continue

        df = pd.read_csv(file_path, sep=None, engine='python')
        df.columns = [col.strip().replace(" ", "") for col in df.columns]

        if 'Protein' not in df.columns:
            continue

        # Use matches_gene_id() to match any protein isoform
        mask = df['Protein'].str.split(';').apply(
            lambda isoforms: any(matches_gene_id(p.strip(), protein_name) for p in isoforms)
        )

        df_protein = df[mask].copy()

        if df_protein.empty:
            sample_cols = [col for col in df.columns if col != 'Protein']
            zero_row = {col: 0 for col in sample_cols}
            zero_row['Protein'] = protein_name
            df_protein = pd.DataFrame([zero_row])

        df_melted = df_protein.melt(id_vars=['Protein'], var_name='Sample', value_name='Intensity')
        df_melted['Band'] = band
        all_data.append(df_melted)

    if not all_data:
        raise ValueError("No data collected. Check protein name or input files.")

    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df


"""Return a list of bands a molecular weight falls into."""
def assign_band(mw_kda):
    return [band for band, low, high in GEL_BANDS if low <= mw_kda <= high]


"""Calculate molecular weight of peptide sequence including water mass."""
def calculate_mw(sequence):
    return (sum(AA_WEIGHTS[aa] for aa in sequence) + WATER_MASS)/1000


def parse_fasta(fasta_path, gene_name):
    sequences = {}
    with open(fasta_path) as f:
        current_id, current_seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and matches_gene_id(current_id, gene_name):
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id and matches_gene_id(current_id, gene_name):
            sequences[current_id] = ''.join(current_seq)
    return sequences

def find_peptide_to_isoform_map(peptide_csv_dir, gene_name):
    isoform_to_peptides = {}
    pattern = fr'^{re.escape(gene_name)}\.'  # Exact match: gene_name followed by dot

    for csv_file in Path(peptide_csv_dir).glob("*.csv"):
        try:
            df = pd.read_csv(csv_file, sep=None, engine="python", encoding="utf-8")
            df.columns = df.columns.str.strip().str.replace('\ufeff', '', regex=False)
            if not {"sequences", "Protein"}.issubset(df.columns):
                continue

            # Filter proteins by exact gene name + dot prefix
            df = df[df['Protein'].str.match(pattern, na=False)]

            for _, row in df.iterrows():
                peptide = row['sequences']
                for isoform in row['Protein'].split(';'):
                    isoform = isoform.strip()
                    if isoform.startswith(f"{gene_name}."):
                        isoform_to_peptides.setdefault(isoform, set()).add(peptide)
        except Exception:
            continue

    return isoform_to_peptides

def get_peptide_exon_mapping(peptides, protein_seq, aa_to_exon):
    mapping = {}
    for pep in peptides:
        match = re.search(re.escape(pep), protein_seq)
        if not match:
            continue
        start, end = match.start(), match.end()
        exons = sorted(set(aa_to_exon.get(i) for i in range(start, end) if i in aa_to_exon))
        mapping[pep] = exons
    return mapping

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def translate_dna(dna_seq):
    return ''.join([GENETIC_CODE.get(dna_seq[i:i+3], 'X') for i in range(0, len(dna_seq) - 2, 3)])
