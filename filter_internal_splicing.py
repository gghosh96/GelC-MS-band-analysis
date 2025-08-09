from constants import BAND_FOLDER
from utils import matches_gene_id
import os
import pandas as pd
from ucsc_api import *
from constants import FASTA_PATH
from utils import translate_dna, parse_fasta, map_exons_to_aa

def filter_proteins_with_min_peptides(protein_list, min_peptides=5):
    """
    Filters the given list of proteins to those having at least `min_peptides` 
    unique peptide sequences across all CSV files in BAND_FOLDER.
    """
    # Dictionary to hold protein -> set of unique peptides
    protein_peptides = {prot: set() for prot in protein_list}

    # Go through all CSV files
    for file in os.listdir(BAND_FOLDER):
        if not file.endswith(".csv"):
            continue
        df = pd.read_csv(os.path.join(BAND_FOLDER, file))

        if 'Protein' not in df.columns or 'sequences' not in df.columns:
            continue
        
        # Iterate through rows
        for _, row in df.iterrows():
            prots = str(row['Protein']).split(';')
            peptide = str(row['sequences']).strip()

            for prot in prots:
                for target_protein in protein_list:
                    if matches_gene_id(prot.strip(), target_protein):
                        protein_peptides[target_protein].add(peptide)

    # Filter to only proteins with ≥ min_peptides unique peptides
    filtered_proteins = [
        prot for prot, peptides in protein_peptides.items()
        if len(peptides) >= min_peptides
    ]

    return filtered_proteins

def flag_splicing_events(canonical_seq, isoform_seq, terminal_threshold=10, internal_length_threshold=10):
    min_len = min(len(canonical_seq), len(isoform_seq))
    first_diff_pos = None
    for i in range(min_len):
        if canonical_seq[i] != isoform_seq[i]:
            first_diff_pos = i
            break

    if first_diff_pos is None:
        if len(canonical_seq) != len(isoform_seq):
            return {
                'terminal_change': True,
                'internal_splicing': False,
                'first_diff_pos': min_len,
                'notes': 'Identical Sequence (Longest Canonical Protein)'
            }
        else:
            return {
                'terminal_change': False,
                'internal_splicing': False,
                'first_diff_pos': None,
                'notes': 'Sequences identical'
            }

    if first_diff_pos < terminal_threshold:
        return {
            'terminal_change': True,
            'internal_splicing': False,
            'first_diff_pos': first_diff_pos,
            'notes': f'First mismatch within first {terminal_threshold} amino acids, likely terminal change'
        }

    suffix_len = 0
    for i in range(1, min_len - first_diff_pos + 1):
        if canonical_seq[-i] != isoform_seq[-i]:
            break
        suffix_len += 1

    canonical_internal_len = len(canonical_seq) - first_diff_pos - suffix_len
    isoform_internal_len = len(isoform_seq) - first_diff_pos - suffix_len

    if canonical_internal_len >= internal_length_threshold or isoform_internal_len >= internal_length_threshold:
        return {
            'terminal_change': False,
            'internal_splicing': True,
            'first_diff_pos': first_diff_pos,
            'notes': f'Internal region length difference >= {internal_length_threshold}, flagged internal splicing'
        }

    return {
        'terminal_change': False,
        'internal_splicing': False,
        'first_diff_pos': first_diff_pos,
        'notes': 'Difference detected but not flagged as terminal or internal splicing'
    }

def flag_isoforms_with_splicing_events(protein_list, terminal_threshold=10, internal_length_threshold=10):
    flagged = {}

    for gene in protein_list:
        print(f"Processing {gene}...")

        # Get canonical transcript info
        tx = fetch_transcript_info(gene)
        if not tx:
            print(f"  No transcript info found for {gene}. Skipping.")
            continue

        # Use map_exons_to_aa to get correct spliced canonical protein
        canonical_protein, aa_ranges, aa_to_exon = map_exons_to_aa(tx)

        print(f"Canonical protein sequence for {gene} (length {len(canonical_protein)} aa):")
        print(canonical_protein)

        # Parse all isoforms from your FASTA file
        isoforms = parse_fasta(FASTA_PATH, gene)

        isoform_flags = []
        for isoform_id, iso_seq in isoforms.items():
            # Include all isoforms without skipping canonical

            result = flag_splicing_events(
                canonical_protein,
                iso_seq,
                terminal_threshold=terminal_threshold,
                internal_length_threshold=internal_length_threshold
            )

            isoform_flags.append({
                'isoform_id': isoform_id,
                'flags': result
            })

        if isoform_flags:
            flagged[gene] = isoform_flags

    return flagged

if __name__ == "__main__":
    protein_list = ["HNRNPA2B1", "MBNL1"]

    # Filter proteins with at least 5 unique peptides
    filtered_proteins = filter_proteins_with_min_peptides(protein_list, min_peptides=5)

    # Flag isoforms with terminal or internal splicing events (and identical ones too)
    flagged_splicing = flag_isoforms_with_splicing_events(filtered_proteins,
                                                          terminal_threshold=10,
                                                          internal_length_threshold=10)

    print("\nProteins with flagged splicing isoforms:")
    for gene, isoforms in flagged_splicing.items():
        print(f"{gene}:")
        for iso in isoforms:
            print(f"  {iso['isoform_id']}: {iso['flags']['notes']}")