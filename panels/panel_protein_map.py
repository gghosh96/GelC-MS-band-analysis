import matplotlib.pyplot as plt

def plot_protein(protein, aa_ranges, gene_name, output_dir, wrap_width=100):
    """
    Visualizes the full protein sequence with exons highlighted.

    Parameters:
    - protein (str): Full amino acid sequence of the protein.
    - aa_ranges (list of tuples): List of (start, end) positions for exons in amino acid coordinates.
    - gene_name (str): Gene symbol to label the plot.
    - output_dir (Path): Directory to save the output image.
    - wrap_width (int): Number of amino acids per line in the wrapped sequence.
    """

    # Calculate how many lines the protein will wrap across based on wrap_width
    num_lines = (len(protein) + wrap_width - 1) // wrap_width

    # Adjust figure height based on the number of wrapped lines
    fig_height = 1.5 + num_lines * 1.2
    fig, ax = plt.subplots(figsize=(max(12, wrap_width * 0.15), fig_height))

    # Define a color palette to differentiate exons
    colors = plt.cm.tab20.colors

    y_base = num_lines  # Starting vertical position for drawing

    legend_seen = set()  # Keep track of exons already added to legend

    # Iterate over each wrapped line of amino acids
    for line_idx in range(num_lines):
        start_idx = line_idx * wrap_width
        end_idx = min((line_idx + 1) * wrap_width, len(protein))
        line_protein = protein[start_idx:end_idx]

        # For each exon, check if it overlaps with the current line and plot its span
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

        # Plot each amino acid character in the line and color it by its exon
        for i, aa in enumerate(line_protein):
            global_pos = start_idx + i
            exon_idx = next((j for j, (start, end) in enumerate(aa_ranges) if start <= global_pos < end), None)
            color = colors[exon_idx % len(colors)] if exon_idx is not None else 'black'
            ax.text(i + 0.5, y_base - line_idx + 0.15, aa,
                    color=color, fontsize=11, fontweight='bold',
                    fontfamily='monospace', ha='center', va='bottom')

    # Set axis and title information representing amino acid coordinates
    ax.set_xlim(0, wrap_width)
    ax.set_ylim(0, y_base + 1)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([])
    ax.set_title(f"{gene_name} — Exon-labeled Protein Sequence")
    ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    plt.tight_layout()

    # Save the final plot to the output directory
    output_file = output_dir / f"{gene_name}_protein_map.png"
    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file
