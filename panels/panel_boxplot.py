import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def plot_panel_a_boxplot(df, protein_name, output_dir):
    # --- Step 1: filter by protein name ---
    if "Gene" in df.columns:
        df = df[df["Gene"] == protein_name]

    if df.empty:
        raise ValueError(f"No data found for protein {protein_name}")

    # --- Step 2: collapse duplicates by summing intensities ---
    group_cols = ["Sample", "Band"]
    if "Gene" in df.columns:
        group_cols.append("Gene")

    df = df.groupby(group_cols, as_index=False)["Intensity"].sum()

    # --- Step 3: assign Tumor/Normal ---
    df['Condition'] = df['Sample'].astype(str).apply(
        lambda x: 'Tumor' if x.endswith('_T') else ('Normal' if x.endswith('_N') else 'Unknown')
    )

    # Keep only Tumor and Normal
    df = df[df['Condition'].isin(['Tumor', 'Normal'])]
    if df.empty:
        raise ValueError("No valid 'Tumor' or 'Normal' samples found in the data.")

    # --- Step 4: normalize intensity ---
    df['RelativeIntensity'] = df['Intensity'] / df['Intensity'].max()
    df = df.dropna(subset=['RelativeIntensity'])

    capped_values = []
    outlier_indices = []

    # --- Step 5: outlier detection & capping per Band × Condition ---
    for (band, condition), group in df.groupby(['Band', 'Condition']):
        vals = group['RelativeIntensity'].copy()
        q1 = np.percentile(vals, 25)
        q3 = np.percentile(vals, 75)
        iqr = q3 - q1
        median = np.median(vals)
        
        # Tukey thresholds
        lower_tukey = q1 - 1.5 * iqr
        upper_tukey = q3 + 1.5 * iqr
        
        # Outlier mask
        outlier_mask = (vals < lower_tukey) | (vals > upper_tukey)
        outlier_indices.extend(group[outlier_mask].index.tolist())
        
        # Cap at median ± 2*IQR
        lower_cap = median - 2 * iqr
        upper_cap = median + 2 * iqr
        vals_capped = vals.clip(lower=lower_cap, upper=upper_cap)
        capped_values.append(vals_capped)
    
    # Replace values
    df['RelativeIntensityCapped'] = pd.concat(capped_values).sort_index()
    df['IsOutlier'] = False
    df.loc[outlier_indices, 'IsOutlier'] = True

    # Split inliers/outliers
    df_inliers = df[~df['IsOutlier']]
    df_outliers = df[df['IsOutlier']]

    # --- Step 6: plotting ---
    sns.set(style="whitegrid")
    plt.figure(figsize=(14, 6))

    hue_order = ['Normal', 'Tumor']

    # Boxplot (capped values, inliers only)
    sns.boxplot(
        data=df_inliers,
        x='Band',
        y='RelativeIntensityCapped',
        hue='Condition',
        hue_order=hue_order,
        palette={'Tumor': '#FF6666', 'Normal': '#87CEFA'},
        dodge=True,
        width=0.6,
        showfliers=False
    )

    # Overlay inlier points
    sns.stripplot(
        data=df_inliers,
        x='Band',
        y='RelativeIntensityCapped',
        hue='Condition',
        hue_order=hue_order,
        palette={'Tumor': '#CC0000', 'Normal': '#4682B4'},
        dodge=True,
        jitter=True,
        marker='o',
        alpha=0.7,
        size=5,
        edgecolor='auto',
        linewidth=0.5
    )

    # Overlay outliers as Xs
    outlier_palette = {'Tumor': "#AF0000", 'Normal': "#0037CF"}
    sns.stripplot(
        data=df_outliers,
        x='Band',
        y='RelativeIntensityCapped',
        hue='Condition',
        hue_order=hue_order,
        palette=outlier_palette,
        dodge=True,
        jitter=True,
        marker='X',
        alpha=0.9,
        size=7,
        edgecolor='white',
        linewidth=1
    )

    # Fix legend
    handles, labels = plt.gca().get_legend_handles_labels()
    n_conditions = len(hue_order)
    plt.legend(handles[:n_conditions], labels[:n_conditions], title='Condition')

    plt.title(f'Relative {protein_name} Intensity Across Gel Bands', fontsize=14)
    plt.xlabel('Band Number')
    plt.ylabel('Relative Intensity')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()

    output_file = output_dir / f"{protein_name}_protein_quant_boxplot.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    return output_file


