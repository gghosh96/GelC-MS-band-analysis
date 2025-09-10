import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Rectangle, Patch
import matplotlib.cm as cm
from utils import *
from ucsc_api import *
from matplotlib.lines import Line2D

def analyze_isoforms_clustered(fasta_path, peptide_csv_dir, gene_name):
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

def collapse_isoform_clusters_with_sas(results, canonical_seq):
    """
    Collapse isoforms into clusters accounting for:
    - Parent isoform (longest / closest to canonical initially)
    - Truncations (may belong to multiple parents)
    - Full-length substitutions (SAS)
    - Truncations of substitutions

    Post-processing:
    1. Redefine cluster parent as heaviest isoform by MW.
    2. If a parent isoform appears as a truncation in another cluster,
       remove its cluster entirely.
    """
    def similarity_to_canonical(seq):
        return sum(a == b for a, b in zip(seq, canonical_seq))

    # --- Detect canonical isoform ---
    canonical_cluster = None
    canonical_parent_id = None
    for iso in results:
        if iso["sequence"] == canonical_seq:
            canonical_parent_id = iso["isoform_id"]
            canonical_cluster = {
                "parent": {**iso, "is_canonical": True},
                "truncations": [],
                "substitutions": [],
                "substituted_truncations": []
            }
            break

    # --- Sort all isoforms by sequence length and similarity ---
    results = sorted(results, key=lambda r: (-len(r["sequence"]), -similarity_to_canonical(r["sequence"])))
    clustered = []
    assigned = set()

    for i, parent in enumerate(results):
        if parent["isoform_id"] in assigned:
            continue
        if canonical_parent_id and parent["isoform_id"] == canonical_parent_id:
            assigned.add(canonical_parent_id)
            clustered.append(canonical_cluster)
            continue

        parent_seq = parent["sequence"]
        cluster = {
            "parent": parent,
            "truncations": [],
            "substitutions": [],
            "substituted_truncations": []
        }

        for j in range(i + 1, len(results)):
            candidate = results[j]
            if candidate["isoform_id"] in assigned:
                continue
            cand_seq = candidate["sequence"]

            if len(cand_seq) == len(parent_seq) and cand_seq != parent_seq:
                cluster["substitutions"].append(candidate)
                assigned.add(candidate["isoform_id"])
            elif cand_seq in parent_seq:
                cluster["truncations"].append(candidate)
            else:
                for sas_var in cluster["substitutions"]:
                    if cand_seq in sas_var["sequence"]:
                        cluster["substituted_truncations"].append(candidate)
                        break

        assigned.add(parent["isoform_id"])

        # --- Reassign heaviest MW as the "true parent" ---
        all_isoforms = [parent] + cluster["substitutions"] + cluster["truncations"] + cluster["substituted_truncations"]
        true_parent = max(all_isoforms, key=lambda x: x["mw_kda"])
        if true_parent["isoform_id"] != parent["isoform_id"]:
            cluster["parent"] = true_parent
            # move old parent into truncations if smaller
            if parent not in cluster["truncations"] and parent["isoform_id"] != true_parent["isoform_id"]:
                cluster["truncations"].append(parent)

        clustered.append(cluster)

    # --- Post-processing: remove clusters whose parent is a truncation elsewhere ---
    truncation_isoforms = set()
    for cluster in clustered:
        for t in cluster["truncations"]:
            truncation_isoforms.add(t["isoform_id"])
        for st in cluster["substituted_truncations"]:
            truncation_isoforms.add(st["isoform_id"])

    filtered_clusters = []
    for c in clustered:
        pid = c["parent"]["isoform_id"]
        if pid in truncation_isoforms and not c["parent"].get("is_canonical", False):
            continue  # remove this cluster
        filtered_clusters.append(c)

    # --- Put canonical cluster at the very start if it exists ---
    if canonical_cluster:
        filtered_clusters = [canonical_cluster] + [c for c in filtered_clusters if c is not canonical_cluster]

    return filtered_clusters

def plot_clustered_isoform_maps_with_sas(results, gene_name, output_dir):
    # Fetch canonical protein and exon mapping
    transcripts = fetch_all_transcripts(gene_name)
    tx = select_longest_cds_transcript(transcripts)
    if tx is None:
        raise ValueError(f"No transcript found for gene {gene_name}")
    canonical_protein, aa_ranges, _ = map_exons_to_aa(tx)
    canonical_protein = canonical_protein.rstrip('*')  # remove stop codons
    aa_ranges = sorted(aa_ranges, key=lambda x: x[0])

    # Collapse isoforms into clusters (post-processing handles true parents)
    clusters = collapse_isoform_clusters_with_sas(results, canonical_protein)

    # Safety check: if no clusters remain, return None
    if not clusters:
        print(f"No clusters remaining for {gene_name} after truncation filtering.")
        return None, None

    # Assign the heaviest MW variant in each cluster as the parent
    for cluster in clusters:
        all_variants = [cluster['parent']] + cluster['truncations'] + cluster['substitutions'] + cluster['substituted_truncations']
        heaviest = max(all_variants, key=lambda x: x['mw_kda'])
        cluster['parent'] = heaviest

    # Prepare exon colors
    colors = cm.get_cmap('tab20').colors
    filtered_colors = [c for c in colors if not (c[0] > 0.8 and c[1] < 0.3 and c[2] < 0.3)]
    exon_color_map = {i: filtered_colors[i % len(filtered_colors)] for i in range(len(aa_ranges))}

    # Try to find canonical sequence in our data
    canonical_cluster = None
    for r in results:
        if r['sequence'] == canonical_protein:
            canonical_cluster = {
                'parent': r,
                'truncations': [],
                'substitutions': [],
                'substituted_truncations': []
            }
            break

    if canonical_cluster:
        clusters.insert(0, canonical_cluster)  # put canonical cluster first

    # Sort remaining clusters by MW for plotting (excluding canonical at top)
    clusters[1:] = sorted(clusters[1:], key=lambda c: c['parent']['mw_kda'])

    # Begin plotting
    fig_height = 0.6 * len(clusters)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    y_offset = 0
    max_len = max(len(c["parent"]["sequence"]) for c in clusters)
    used_exons = set()
    txt_lines = []

    for idx, cluster in enumerate(clusters):
        parent = cluster["parent"]
        iso_seq = parent["sequence"]
        iso_y = y_offset + 1

        # Draw parent line (dotted for all clusters)
        line = Line2D([0, len(iso_seq)], [iso_y, iso_y], color='black', linewidth=4)
        line.set_dashes([1, 0.5])
        ax.add_line(line)

        # Map exons
        for exon_idx, (start, end) in enumerate(aa_ranges):
            exon_seq = canonical_protein[start:end]
            pos_in_iso = iso_seq.find(exon_seq)
            if pos_in_iso != -1:
                ax.hlines(iso_y, pos_in_iso, pos_in_iso + len(exon_seq),
                          colors=exon_color_map[exon_idx], linewidth=4)
                used_exons.add(exon_idx)

        # Overlay peptides
        for pep in parent["peptides"]:
            pep_start = iso_seq.find(pep)
            if pep_start == -1:
                continue
            # find peptide position in canonical (ignore *)
            pep_clean = pep.replace('*', '')
            match = re.search(re.escape(pep_clean), canonical_protein)
            if not match:
                ax.add_patch(Rectangle((pep_start, iso_y-0.2), len(pep), 0.4,
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
            # split peptide by exon boundaries
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
                ax.add_patch(Rectangle((offset_pep, iso_y-0.2), len(seg), 0.4,
                                       facecolor=color, edgecolor='black', lw=0.5, zorder=3))
                offset_pep += len(seg)

        # Label cluster
        n_trunc = len(cluster["truncations"])
        n_sas = len(cluster["substitutions"])
        n_sub_trunc = len(cluster["substituted_truncations"])
        label_name = "Canonical" if idx == 0 and canonical_cluster else parent['isoform_id']
        label = f"{label_name}\n{parent['mw_kda']} kDa, Bands: {', '.join(map(str,parent['bands']))}"
        extras = []
        if n_trunc: extras.append(f"{n_trunc} truncations")
        if n_sas: extras.append(f"{n_sas} substitutions")
        if n_sub_trunc: extras.append(f"{n_sub_trunc} substituted truncations")
        if extras:
            label += "\n(+ " + ", + ".join(extras) + ")"
        ax.text(-0.01*max_len, iso_y, label, va="center", ha="right", fontsize=8)

        # TXT output
        all_isoforms = [parent["isoform_id"]] + \
                       [c["isoform_id"] for c in cluster["truncations"]] + \
                       [c["isoform_id"] for c in cluster["substitutions"]] + \
                       [c["isoform_id"] for c in cluster["substituted_truncations"]]
        txt_lines.append(f"Main isoform: {parent['isoform_id']} -> " +
                         ", ".join([i for i in all_isoforms if i != parent["isoform_id"]]))

        y_offset += 1.5

    ax.set_xlim(0, max_len)
    ax.set_ylim(0, y_offset+1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"Detected Peptides in {gene_name} Isoforms (Clustered + SAS)")

    # Legend
    exon_legend_patches = [Patch(facecolor=exon_color_map[i], label=f"Exon {i+1}") 
                           for i in range(len(aa_ranges)) if i in used_exons]
    ax.legend(handles=exon_legend_patches, title="Exons", bbox_to_anchor=(1.02,1),
              loc='upper left', borderaxespad=0.)

    plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{gene_name}_clustered_protein_isoform_map_SAS.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Write TXT
    txt_file = output_dir / f"{gene_name}_clustered_isoform_map_clusters.txt"
    with open(txt_file, "w") as f:
        f.write("\n".join(txt_lines))

    return output_file, txt_file
