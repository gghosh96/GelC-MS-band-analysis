import matplotlib.pyplot as plt
import seaborn as sns

# Plot relative protein quantification boxplot across gel bands for a given protein
def plot_panel_a_boxplot(df, protein_name, output_dir):
    
    # Assign sample condition based on sample name suffix (_T for Tumor, _N for Normal)
    df['Condition'] = df['Sample'].astype(str).apply(
        lambda x: 'Tumor' if x.endswith('_T') else ('Normal' if x.endswith('_N') else 'Unknown')
    )
    
    # Filter dataframe to keep only Tumor and Normal samples
    df = df[df['Condition'].isin(['Tumor', 'Normal'])]
    
    # Raise error if no valid Tumor or Normal samples found
    if df.empty:
        raise ValueError("No valid 'Tumor' or 'Normal' samples found in the data.")
    
    # Calculate relative intensity by normalizing intensity values within the dataset
    df['RelativeIntensity'] = df['Intensity'] / df['Intensity'].max()
    
    # Drop rows where relative intensity could not be computed (e.g., missing data)
    df = df.dropna(subset=['RelativeIntensity'])

# Create boxplot showing distribution of relative intensities per gel band, separated by condition
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

    # Overlay individual data points to show spread and overlap, colored by condition
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
    output_file = output_dir /  f"{protein_name}_protein_quant_boxplot.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    return output_file
