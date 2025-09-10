from constants import BAND_FOLDER
from utils import matches_gene_id
import os
import pandas as pd
from ucsc_api import *
from constants import FASTA_PATH
from utils import translate_dna, parse_fasta, map_exons_to_aa

def build_protein_peptide_dict(protein_list):
    """
    Collects all unique peptides detected per protein across all CSVs in BAND_FOLDER.
    """
    protein_peptides = {prot: set() for prot in protein_list}

    for file in os.listdir(BAND_FOLDER):
        if not file.endswith(".csv"):
            continue
        df = pd.read_csv(os.path.join(BAND_FOLDER, file))

        if 'Protein' not in df.columns or 'sequences' not in df.columns:
            continue
        
        for _, row in df.iterrows():
            prots = str(row['Protein']).split(';')
            peptide = str(row['sequences']).strip()

            for prot in prots:
                for target_protein in protein_list:
                    if matches_gene_id(prot.strip(), target_protein):
                        protein_peptides[target_protein].add(peptide)

    return protein_peptides

def filter_proteins_by_coverage(protein_list, protein_peptides, coverage_threshold=0.3):
    """
    Filters proteins based on fraction of canonical sequence covered by detected peptides.
    """
    filtered_proteins = []

    for protein in protein_list:
        transcripts = fetch_all_transcripts(protein)
        tx = select_longest_cds_transcript(transcripts)
        if not tx:
            continue

        canonical_protein, _, _ = map_exons_to_aa(tx)
        covered = [0] * len(canonical_protein)

        for peptide in protein_peptides.get(protein, []):
            start = canonical_protein.find(peptide)
            if start != -1:
                for i in range(start, start + len(peptide)):
                    covered[i] = 1

        coverage_fraction = sum(covered) / len(canonical_protein)
        if coverage_fraction >= coverage_threshold:
            filtered_proteins.append(protein)

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

        transcripts = fetch_all_transcripts(gene)
        tx = select_longest_cds_transcript(transcripts)
        if not tx:
            print(f"  No transcript info found for {gene}. Skipping.")
            continue

        canonical_protein, aa_ranges, aa_to_exon = map_exons_to_aa(tx)

        print(f"Canonical protein sequence for {gene} (length {len(canonical_protein)} aa):")
        print(canonical_protein)

        isoforms = parse_fasta(FASTA_PATH, gene)

        isoform_flags = []
        for isoform_id, iso_seq in isoforms.items():
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

    # Build protein -> peptide mapping
    protein_peptides = build_protein_peptide_dict(protein_list)

    # Filter proteins by ≥30% sequence coverage
    filtered_proteins = filter_proteins_by_coverage(protein_list, protein_peptides, coverage_threshold=0.3)

    # Flag isoforms with terminal or internal splicing events
    flagged_splicing = flag_isoforms_with_splicing_events(filtered_proteins,
                                                          terminal_threshold=10,
                                                          internal_length_threshold=10)

    print("\nProteins with flagged splicing isoforms:")
    for gene, isoforms in flagged_splicing.items():
        print(f"{gene}:")
        for iso in isoforms:
            print(f"  {iso['isoform_id']}: {iso['flags']['notes']}")
