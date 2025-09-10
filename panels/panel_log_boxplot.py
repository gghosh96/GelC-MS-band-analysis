import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def plot_log_protein_boxplot(df, protein_name, output_dir):
    """
    Create a log10-scaled protein quantification boxplot with Tukey-style outlier detection.
    Highlights subtle differences in lower molecular weight bands,
    using 10^3 as baseline, and performs Tumor vs Normal t-tests.
    """

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

    # --- Step 4: apply baseline (10^3) and log transform ---
    df['IntensityAdjusted'] = df['Intensity'].clip(lower=1e3)
    df['LogIntensity'] = np.log10(df['IntensityAdjusted'])

    # --- Step 5: run t-tests per band (Tumor vs Normal) ---
    ttest_results = []
    for band, group in df.groupby("Band"):
        tumor_vals = group.loc[group['Condition'] == 'Tumor', 'LogIntensity']
        normal_vals = group.loc[group['Condition'] == 'Normal', 'LogIntensity']
        if len(tumor_vals) > 1 and len(normal_vals) > 1:
            t_stat, p_val = ttest_ind(tumor_vals, normal_vals, equal_var=False)
            ttest_results.append({
                "Protein": protein_name,
                "Band": band,
                "T-stat": t_stat,
                "P-value": p_val
            })

    ttest_df = pd.DataFrame(ttest_results)
    output_stats_file = None
    if not ttest_df.empty:
        output_stats_file = output_dir / f"{protein_name}_log_ttest_results.csv"
        ttest_df.to_csv(output_stats_file, index=False)

    # --- Step 6: Tukey outlier detection and capping per Band × Condition ---
    capped_values = []
    outlier_indices = []

    for (band, condition), group in df.groupby(['Band', 'Condition']):
        vals = group['LogIntensity'].copy()
        q1 = np.percentile(vals, 25)
        q3 = np.percentile(vals, 75)
        iqr = q3 - q1
        median = np.median(vals)
        
        lower_tukey = q1 - 1.5 * iqr
        upper_tukey = q3 + 1.5 * iqr
        
        outlier_mask = (vals < lower_tukey) | (vals > upper_tukey)
        outlier_indices.extend(group[outlier_mask].index.tolist())
        
        # Cap inliers to median ± 2*IQR
        lower_cap = median - 2 * iqr
        upper_cap = median + 2 * iqr
        vals_capped = vals.clip(lower=lower_cap, upper=upper_cap)
        capped_values.append(vals_capped)

    df['LogIntensityCapped'] = pd.concat(capped_values).sort_index()
    df['IsOutlier'] = False
    df.loc[outlier_indices, 'IsOutlier'] = True

    df_inliers = df[~df['IsOutlier']]
    df_outliers = df[df['IsOutlier']]

    # --- Step 7: plotting ---
    sns.set(style="whitegrid")
    plt.figure(figsize=(14, 6))

    hue_order = ['Normal', 'Tumor']

    # Boxplot (inliers only)
    sns.boxplot(
        data=df_inliers,
        x='Band',
        y='LogIntensityCapped',
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
        y='LogIntensityCapped',
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
        y='LogIntensityCapped',
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

    plt.title(f'Log10 {protein_name} Intensity Across Gel Bands', fontsize=14)
    plt.xlabel('Band Number')
    plt.ylabel('Log10(Intensity), baseline=1e3')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()

    output_file = output_dir / f"{protein_name}_log_protein_quant_boxplot.png"
    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file, output_stats_file

