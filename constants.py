# === Fixed Configuration ===
FOLDER_PATH = "Gel_Band_Protein_Quants/"
BAND_RANGE = range(1, 17)
SAVE_PLOTS = True
BAND_FOLDER = "Gel_Band_Peptide_Quants/"
FASTA_PATH = "Transcriptome_peptides_final_headers.fasta"
PEPTIDE_CSV_DIR = BAND_FOLDER

# === Genetic Code Dictionary ===
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

# === Gel Band MW Ranges (kDa) ===
GEL_BANDS = [
    (1, 8, 14), (2, 13, 17), (3, 15, 19), (4, 18, 22), (5, 20, 25),
    (6, 24, 30), (7, 28, 35), (8, 32, 40), (9, 35, 44), (10, 42, 55),
    (11, 47, 63), (12, 58, 90), (13, 75, 160), (14, 140, 280),
    (15, 240, 450), (16, 450, float("inf"))
]

# === Amino Acid Molecular Weights (Da) ===
AA_WEIGHTS = {
    'A': 71.0788,  'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
    'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
    'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
    'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
    'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
}

WATER_MASS = 18.01528