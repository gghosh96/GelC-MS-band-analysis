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
from PIL import Image

# === Fixed Configuration ===
FOLDER_PATH = "Gel_Band_Protein_Quants/"
BAND_RANGE = range(1, 17)
SAVE_PLOTS = True
BAND_FOLDER = "Gel_Band_Peptide_Quants/"
FASTA_PATH = "Transcriptome_peptides_final_headers.fasta"
PEPTIDE_CSV_DIR = BAND_FOLDER

# Genetic code dictionary
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}
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


def plot_panel_a_boxplot(df, protein_name):
    df['Condition'] = df['Sample'].astype(str).apply(
        lambda x: 'Tumor' if x.endswith('_T') else ('Normal' if x.endswith('_N') else 'Unknown')
    )
    df = df[df['Condition'].isin(['Tumor', 'Normal'])]
    if df.empty:
        raise ValueError("No valid 'Tumor' or 'Normal' samples found in the data.")
    df['RelativeIntensity'] = df['Intensity'] / df['Intensity'].max()
    df = df.dropna(subset=['RelativeIntensity'])
    sns.set(style="whitegrid")
    plt.figure(figsize=(14, 6))
    sns.boxplot(
        data=df,
        x='Band',
        y='RelativeIntensity',
        hue='Condition',
        palette={'Tumor': '#FF6666', 'Normal': '#87CEFA'},
        dodge=True,
        width=0.6,
        showfliers=False
    )
    sns.stripplot(
        data=df,
        x='Band',
        y='RelativeIntensity',
        hue='Condition',
        palette={'Tumor': '#CC0000', 'Normal': '#4682B4'},
        dodge=True,
        jitter=True,
        marker='o',
        alpha=0.7,
        size=5,
        edgecolor='auto',
        linewidth=0.5
    )
    handles, labels = plt.gca().get_legend_handles_labels()
    n = int(len(handles) / 2)
    plt.legend(handles[:n], labels[:n], title='Condition')
    plt.title(f'Relative {protein_name} Intensity Across Gel Bands', fontsize=14)
    plt.xlabel('Band Number')
    plt.ylabel('Relative Intensity (0–1)')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    output_file = f"{protein_name}_panelA_boxplot.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    return output_file

def connect_ucsc():
    return mysql.connector.connect(user='genome', host='genome-mysql.soe.ucsc.edu', database='hg38')

def fetch_transcript_info(gene_name):
    conn = connect_ucsc()
    cur = conn.cursor(dictionary=True)
    query = """
        SELECT * FROM refGene
        WHERE name2 = %s AND cdsStart < cdsEnd
        ORDER BY exonCount DESC LIMIT 1
    """
    cur.execute(query, (gene_name,))
    row = cur.fetchone()
    conn.close()
    return row
    
def fetch_genomic_sequence(chrom, start, end):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    res = requests.get(url)
    res.raise_for_status()
    return res.json()["dna"].upper()

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def translate_dna(dna_seq):
    return ''.join([GENETIC_CODE.get(dna_seq[i:i+3], 'X') for i in range(0, len(dna_seq) - 2, 3)])

def map_exons_to_aa(transcript):
    exon_starts = list(map(int, transcript['exonStarts'].decode().strip(',').split(',')))
    exon_ends = list(map(int, transcript['exonEnds'].decode().strip(',').split(',')))
    cds_start, cds_end = transcript['cdsStart'], transcript['cdsEnd']
    chrom, strand = transcript['chrom'], transcript['strand']

    coding_dna = ""
    aa_ranges, current_aa_pos = [], 0

    for s, e in zip(exon_starts, exon_ends):
        exon_start = max(s, cds_start)
        exon_end = min(e, cds_end)
        if exon_start >= exon_end:
            continue
        seq = fetch_genomic_sequence(chrom, exon_start, exon_end)
        coding_dna += seq
        aa_len = (exon_end - exon_start) // 3
        aa_ranges.append((current_aa_pos, current_aa_pos + aa_len))
        current_aa_pos += aa_len

    if strand == '-':
        coding_dna = reverse_complement(coding_dna)
        total_len = current_aa_pos
        aa_ranges = [(total_len - end, total_len - start) for (start, end) in reversed(aa_ranges)]

    protein = translate_dna(coding_dna)
    aa_to_exon = {i: exon_idx for exon_idx, (start, end) in enumerate(aa_ranges) for i in range(start, end)}

    return protein, aa_ranges, aa_to_exon

def plot_protein(protein, aa_ranges, gene_name, wrap_width=100):
    num_lines = (len(protein) + wrap_width - 1) // wrap_width
    fig_height = 1.5 + num_lines * 1.2
    fig, ax = plt.subplots(figsize=(max(12, wrap_width * 0.15), fig_height))
    colors = plt.cm.tab20.colors
    y_base = num_lines
    legend_seen = set()
    for line_idx in range(num_lines):
        start_idx = line_idx * wrap_width
        end_idx = min((line_idx + 1) * wrap_width, len(protein))
        line_protein = protein[start_idx:end_idx]
        for i, (exon_start, exon_end) in enumerate(aa_ranges):
            seg_start = max(start_idx, exon_start)
            seg_end = min(end_idx, exon_end)
            if seg_start < seg_end:
                x0 = seg_start - start_idx
                x1 = seg_end - start_idx
                label = f"Exon {i+1}" if i not in legend_seen else None
                if label:
                    legend_seen.add(i)
                ax.plot([x0 + 0.5, x1 + 0.5], [y_base - line_idx] * 2, lw=10,
                        color=colors[i % len(colors)], label=label)
        for i, aa in enumerate(line_protein):
            global_pos = start_idx + i
            exon_idx = next((j for j, (start, end) in enumerate(aa_ranges) if start <= global_pos < end), None)
            color = colors[exon_idx % len(colors)] if exon_idx is not None else 'black'
            ax.text(i + 0.5, y_base - line_idx + 0.15, aa,
                    color=color, fontsize=11, fontweight='bold',
                    fontfamily='monospace', ha='center', va='bottom')
    ax.set_xlim(0, wrap_width)
    ax.set_ylim(0, y_base + 1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"{gene_name} — Exon-labeled Protein Sequence")
    ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    plt.tight_layout()
    output_file = f"{gene_name}_exon_colored_protein.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    return output_file

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

def plot_panel_b(protein_name):
    band_files = sorted(
        [(int(re.search(r'Band_(\d+)', f).group(1)), f) for f in os.listdir(BAND_FOLDER) if f.endswith('.csv')],
        key=lambda x: x[0]
    )
    bands = [str(band_num) for band_num, _ in band_files]
    peptide_data = {}

    # Collect peptide data with exact matching
    for band_num, filename in band_files:
        path = os.path.join(BAND_FOLDER, filename)
        df = pd.read_csv(path)
        df['Protein'] = df['Protein'].astype(str)

        # Keep rows where any semicolon-separated isoform starts with protein_name + "."
        df = df[df['Protein'].str.split(';').apply(lambda prots: any(matches_gene_id(p.strip(), protein_name) for p in prots))]

        if df.empty:
            continue

        intensity_cols = [col for col in df.columns if '_T' in col or '_N' in col]
        max_intensity = df[intensity_cols].max(axis=1).replace(0, np.nan)

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

    tx = fetch_transcript_info(protein_name)
    protein_seq, aa_ranges, aa_to_exon = map_exons_to_aa(tx)
    peptides = sorted(peptide_data, key=lambda p: protein_seq.find(p) if protein_seq.find(p) != -1 else float('inf'))
    pep_to_exons = get_peptide_exon_mapping(peptides, protein_seq, aa_to_exon)
    exon_colors = plt.cm.tab20.colors

    n_rows = len(peptides)
    n_cols = len(bands)

    fig_width = max(12, n_cols * 1.3)
    fig_height = max(4, n_rows * 0.7)
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(n_rows, n_cols, figure=fig, wspace=0.05, hspace=0.05)

    all_fc = [data[band]['fold_change'] for data in peptide_data.values() for band in data]
    max_abs = max(abs(x) for x in all_fc) if all_fc else 1

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
        segs, current_exon, current_seg = [], exon_idxs[0], pep[0]
        for aa_idx, exon in zip(range(1, len(pep)), exon_idxs[1:]):
            if exon == current_exon:
                current_seg += pep[aa_idx]
            else:
                segs.append((current_seg, current_exon))
                current_seg, current_exon = pep[aa_idx], exon
        segs.append((current_seg, current_exon))
        offset_x = label_x
        for seg, exon in segs:
            color = exon_colors[exon % len(exon_colors)] if exon is not None else 'gray'
            fig.text(offset_x, y_pos, seg, va='center', ha='left',
                     fontsize=16, fontweight='bold', fontfamily='monospace', color=color)
            offset_x += len(seg) * char_width + gap

    for j, band in enumerate(bands):
        fig.text((j + 0.5) / n_cols, -0.01, f'Band {band}', va='top', ha='center', fontsize=12)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

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

    title_offset = (fig_height + 0.4) / fig_height
    fig.suptitle(f"{protein_name} Peptide Heatmap (Tumor vs Normal)",
                 fontsize=18, y=title_offset)

    output_file = f"{protein_name}_panelB_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_file

# === Constants ===
GEL_BANDS = [
    (1, 8, 14), (2, 13, 17), (3, 15, 19), (4, 18, 22), (5, 20, 25),
    (6, 24, 30), (7, 28, 35), (8, 32, 40), (9, 35, 44), (10, 42, 55),
    (11, 47, 63), (12, 58, 90), (13, 75, 160), (14, 140, 280),
    (15, 240, 450), (16, 450, float("inf"))
]

def assign_band(mw_kda):
    return [band for band, low, high in GEL_BANDS if low <= mw_kda <= high]

def calculate_mw(sequence):
    aa_weights = {
        'A': 71.0788,  'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
        'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
        'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
        'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
        'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }
    water_mass = 18.01528
    return (sum(aa_weights.get(aa, 0.0) for aa in sequence) + water_mass) / 1000

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

def analyze_isoforms(fasta_path, peptide_csv_dir, gene_name):
    fasta_seqs = parse_fasta(fasta_path, gene_name)
    isoform_to_peptides = find_peptide_to_isoform_map(peptide_csv_dir, gene_name)
    results = []

    for isoform_id in isoform_to_peptides:
        matching_keys = [k for k in fasta_seqs if isoform_id in k]
        for key in matching_keys:
            seq = fasta_seqs[key]
            mw = calculate_mw(seq)
            bands = assign_band(mw)
            peptides = sorted(isoform_to_peptides[isoform_id])
            results.append({
                "isoform_id": key,
                "mw_kda": round(mw, 2),
                "bands": bands,
                "sequence": seq,
                "peptides": peptides
            })

    results.sort(key=lambda r: r["mw_kda"])
    return results

# Note: fetch_transcript_info and map_exons_to_aa are assumed to be defined elsewhere

def plot_isoform_maps(results, gene_name):
    tx = fetch_transcript_info(gene_name)
    canonical_protein, _, aa_to_exon = map_exons_to_aa(tx)
    exon_colors = plt.cm.tab20.colors

    fig_height = 0.6 * len(results)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    y_offset = 0
    max_len = max(len(r["sequence"]) for r in results)
    used_exons = set()

    for r in results:
        seq_len = len(r["sequence"])
        iso_y = y_offset + 1
        ax.hlines(iso_y, 0, seq_len, color="black", linewidth=2)

        for pep in r["peptides"]:
            pep_start = r["sequence"].find(pep)
            if pep_start == -1:
                continue
            match = re.search(re.escape(pep), canonical_protein)
            if not match:
                ax.add_patch(plt.Rectangle((pep_start, iso_y - 0.2), len(pep), 0.4,
                                           facecolor='gray', edgecolor='black', lw=0.5))
                continue

            start, end = match.start(), match.end()
            exon_indices = [aa_to_exon.get(pos, None) for pos in range(start, end)]

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

            offset = pep_start
            for seg, exon in label_segments:
                color = exon_colors[exon % len(exon_colors)] if exon is not None else 'gray'
                if exon is not None:
                    used_exons.add(exon)
                ax.add_patch(plt.Rectangle((offset, iso_y - 0.2), len(seg), 0.4,
                                           facecolor=color, edgecolor='black', lw=0.5))
                offset += len(seg)

        label = f"{r['isoform_id']}\n{r['mw_kda']} kDa, Bands: {', '.join(map(str, r['bands']))}"
        ax.text(-0.01 * max_len, iso_y, label, va="center", ha="right", fontsize=8)
        y_offset += 1.5

    ax.set_xlim(0, max_len)
    ax.set_ylim(0, y_offset + 1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"Detected Peptides in {gene_name} Isoforms (Exon-labeled)")

    # === Add Exon Legend ===
    from matplotlib.patches import Patch
    exon_legend_patches = [
        Patch(facecolor=exon_colors[exon % len(exon_colors)], label=f"Exon {exon}")
        for exon in sorted(used_exons)
    ]
    ax.legend(handles=exon_legend_patches, title="Exons", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

    plt.tight_layout()
    output_file = f"{gene_name}_panelD_isoform_map.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_file


def analyze_all_panels(gene_name):
    print(f"Loading data and plotting Panel A for {gene_name}...")
    df = load_diffacto_data(gene_name)
    panelA_file = plot_panel_a_boxplot(df, gene_name)

    print(f"Fetching transcript and plotting exon-colored protein for {gene_name}...")
    tx = fetch_transcript_info(gene_name)
    protein_seq, aa_ranges, aa_to_exon = map_exons_to_aa(tx)
    panelProtein_file = plot_protein(protein_seq, aa_ranges, gene_name)

    print(f"Plotting Panel B heatmap for {gene_name}...")
    panelB_file = plot_panel_b(gene_name)

    print(f"Analyzing isoforms and plotting isoform peptide maps for {gene_name}...")
    isoform_results = analyze_isoforms(FASTA_PATH, PEPTIDE_CSV_DIR, gene_name)
    isoform_file = plot_isoform_maps(isoform_results, gene_name)

    print("All panels generated and saved individually.")
    
    # Return a dict or list of individual files
    return {
        'panelA': panelA_file,
        'panelProtein': panelProtein_file,
        'panelB': panelB_file,
        'isoform': isoform_file
    }

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    from PIL import Image

    if len(sys.argv) < 2:
        print("Usage: python protein_panel_analysis.py <GENE_NAME>")
        sys.exit(1)

    gene = sys.argv[1]
    panel_files = analyze_all_panels(gene)

    for panel_name, filepath in panel_files.items():
        print(f"Displaying {panel_name} saved at: {filepath}")
        img = Image.open(filepath)
        img.show()