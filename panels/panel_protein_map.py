import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_protein(protein, aa_ranges, gene_name, output_dir, wrap_width=100):
    import matplotlib.cm as cm
    import matplotlib.patches as patches

    colors = cm.get_cmap('tab20').colors

    def is_red(color, threshold=0.8):
        r, g, b = color[:3]
        return (r > threshold) and (g < 0.3) and (b < 0.3)

    filtered_colors = [c for c in colors if not is_red(c)]

    # Ensure exons are in N-terminal → C-terminal order
    aa_ranges = sorted(aa_ranges, key=lambda x: x[0])

    # Validate exon coordinates
    for i, (s, e) in enumerate(aa_ranges):
        assert s < e, f"Exon {i+1} has start >= end: {(s, e)}"

    num_lines = (len(protein) + wrap_width - 1) // wrap_width
    fig_height = 1.5 + num_lines * 1.2
    fig_width = max(12, wrap_width * 0.15)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    y_base = num_lines

    # Predefine legend handles in exon order
    legend_handles = [patches.Patch(color=filtered_colors[i % len(filtered_colors)],
                                    label=f"Exon {i+1}") for i in range(len(aa_ranges))]

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
                color = filtered_colors[i % len(filtered_colors)]
                ax.plot([x0 + 0.5, x1 + 0.5], [y_base - line_idx] * 2,
                        lw=10, color=color)

        for i, aa in enumerate(line_protein):
            global_pos = start_idx + i
            exon_idx = next((j for j, (s, e) in enumerate(aa_ranges) if s <= global_pos < e), None)
            color = filtered_colors[exon_idx % len(filtered_colors)] if exon_idx is not None else 'black'
            ax.text(i + 0.5, y_base - line_idx + 0.15, aa,
                    color=color, fontsize=11, fontweight='bold',
                    fontfamily='monospace', ha='center', va='bottom')

    ax.set_xlim(0, wrap_width)
    ax.set_ylim(0, y_base + 1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"{gene_name} — Exon-labeled Protein Sequence")
    ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1.15, 1))
    plt.tight_layout()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{gene_name}_protein_map.png"
    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file
