import sys
import os
from protein_panel_analysis import analyze_all_panels
from utils import matches_gene_id, parse_fasta, translate_dna, map_exons_to_aa
from filter_internal_splicing import (
    flag_splicing_events,
    fetch_genomic_sequence,
    build_protein_peptide_dict,      # Add this
    filter_proteins_by_coverage       # Add this
)
from ucsc_api import fetch_all_transcripts, select_longest_cds_transcript

TERMINAL_THRESHOLD = 10
INTERNAL_LENGTH_THRESHOLD = 10
COVERAGE_THRESHOLD = 0.3  # 30% coverage

def analyze_multiple_genes(gene_list):
    results = {}

    # Build peptide dictionary for all genes upfront
    peptide_dict = build_protein_peptide_dict(gene_list)

    for gene in gene_list:
        print(f"\n=== Analyzing {gene} ===")  # Always print this first

        # Filter based on 30% coverage instead of 5 peptides
        filtered_genes = filter_proteins_by_coverage([gene], peptide_dict, coverage_threshold=COVERAGE_THRESHOLD)
        if not filtered_genes:
            print(f"Skipping {gene}: less than {COVERAGE_THRESHOLD*100:.0f}% coverage.")
            continue

        gene_result = {}

        # Panel analysis
        try:
            panel_files = analyze_all_panels(gene)
            gene_result['panels'] = panel_files
        except Exception as e:
            print(f"Error analyzing panels for {gene}: {e}")
            gene_result['panels'] = None

        # Splicing event flagging
        try:
            transcripts = fetch_all_transcripts(gene)
            tx = select_longest_cds_transcript(transcripts)
            if not tx:
                print(f"  No transcript info found for {gene}. Skipping splicing analysis.")
                gene_result['splicing_flags'] = None
            else:
                canonical_protein, aa_ranges, aa_to_exon = map_exons_to_aa(tx)

                isoforms = parse_fasta('Transcriptome_peptides_final_headers.fasta', gene)
                isoform_flags = []
                for isoform_id, iso_seq in isoforms.items():
                    if isoform_id == f"{gene}.1":
                        continue  # skip canonical isoform

                    flags = flag_splicing_events(
                        canonical_protein,
                        iso_seq,
                        terminal_threshold=TERMINAL_THRESHOLD,
                        internal_length_threshold=INTERNAL_LENGTH_THRESHOLD,
                    )
                    isoform_flags.append({
                        'isoform_id': isoform_id,
                        'flags': flags
                    })

                gene_result['splicing_flags'] = isoform_flags if isoform_flags else None

                # Prepare output text
                output_lines = []
                output_lines.append(f"Canonical protein sequence for {gene} (length {len(canonical_protein)} aa):\n")
                output_lines.append(canonical_protein + "\n\n")

                if isoform_flags:
                    output_lines.append("Splicing events flagged for isoforms:\n")
                    for iso in isoform_flags:
                        iso_id = iso['isoform_id']
                        notes = iso['flags']['notes']
                        output_lines.append(f"Isoform {iso_id}: {notes}\n")
                else:
                    output_lines.append("No splicing events flagged.\n")

                output_dir = os.path.join("data", gene)
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, f"{gene}_isoform_classification.txt")
                with open(output_path, "w") as f:
                    f.writelines(output_lines)

                print("Isoform classification txt file created!")

        except Exception as e:
            print(f"Error flagging splicing events for {gene}: {e}")
            gene_result['splicing_flags'] = None

        results[gene] = gene_result

    return results

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_batch.py <GENE_NAME_1> <GENE_NAME_2> ...")
        sys.exit(1)

    genes = sys.argv[1:]
    analyze_multiple_genes(genes)