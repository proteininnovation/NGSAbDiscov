# IPI 
# Antibody Developability Prediction Platform and Utilities service
# Sequence clustering and heatmap visualization utility
# Design : Hoan Nguyen

## pip install pandas Levenshtein numpy scipy biopython matplotlib seaborn


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from Levenshtein import distance as levenshtein_distance
from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# ====================== FUNCTIONS ======================

def compute_distance_matrix(sequences, method='levenshtein'):
    n = len(sequences)
    dist_matrix = np.zeros((n, n))

    if method == 'levenshtein':
        for i in range(n):
            for j in range(i + 1, n):
                d = levenshtein_distance(sequences[i], sequences[j])
                dist_matrix[i, j] = d
                dist_matrix[j, i] = d

    elif method == 'blosum62':
        max_score = 11
        for i in range(n):
            seq1 = sequences[i]
            len1 = len(seq1)
            for j in range(i + 1, n):
                seq2 = sequences[j]
                len2 = len(seq2)
                alignment_length = max(len1, len2)
                score = 0
                padded1 = seq1.ljust(alignment_length)
                padded2 = seq2.ljust(alignment_length)
                for a, b in zip(padded1, padded2):
                    score += BLOSUM62.get((a, b), BLOSUM62.get((b, a), -4))
                normalized_score = score / (alignment_length * max_score)
                dist = max(0, 1 - normalized_score)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

    return dist_matrix

def create_custom_label(row, label_columns, separator=':'):
    """Combine multiple columns into one label, handling missing values."""
    parts = []
    for col in label_columns:
        val = row.get(col, '')
        if pd.isna(val) or str(val).strip() in ['', 'nan', 'NaN']:
            val = 'NA'
        else:
            val = str(val).strip()
        parts.append(val)
    return separator.join(parts)

def plot_dendrogram_and_heatmap(df, sequence_col, title_prefix="",
                                distance_method='levenshtein', color_threshold=None,
                                label_columns=['BARCODE'], label_separator=':',
                                save_prefix=None):
    if len(df) < 2:
        print(f"Skipping {title_prefix}: not enough sequences.")
        return

    sequences = df[sequence_col].astype(str).tolist()
    
    # Create custom combined labels
    df['custom_label'] = df.apply(
        lambda row: create_custom_label(row, label_columns, label_separator), axis=1
    )
    labels = df['custom_label'].tolist()

    print(f"Computing {distance_method} distances for {title_prefix} ({len(df)} sequences)...")
    dist_matrix = compute_distance_matrix(sequences, method=distance_method)
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')

    order = leaves_list(Z)
    dist_ordered = dist_matrix[np.ix_(order, order)]
    labels_ordered = [labels[i] for i in order]

    # Dendrogram
    fig_width = max(8, len(df) * 0.5)
    plt.figure(figsize=(fig_width, 6))
    dendrogram(Z, labels=labels_ordered, leaf_rotation=90, leaf_font_size=9,
               color_threshold=color_threshold, above_threshold_color='grey')
    plt.title(f"{title_prefix} Dendrogram\n({sequence_col} - {distance_method.capitalize()})", fontsize=14)
    plt.ylabel("Distance")
    plt.tight_layout()
    if save_prefix:
        plt.savefig(f"{save_prefix}_dendrogram.png", dpi=300, bbox_inches='tight')
        print(f"Saved: {save_prefix}_dendrogram.png")
    plt.show()

    # Heatmap
    plt.figure(figsize=(8 + len(df)*0.15, 8))
    sns.heatmap(dist_ordered, cmap="YlOrRd_r", xticklabels=labels_ordered, yticklabels=labels_ordered,
                cbar_kws={"label": "Distance"})
    plt.title(f"{title_prefix} Distance Heatmap\n({sequence_col} - {distance_method.capitalize()})")
    plt.xticks(rotation=90)
    plt.tight_layout()
    if save_prefix:
        plt.savefig(f"{save_prefix}_heatmap.png", dpi=300, bbox_inches='tight')
        print(f"Saved: {save_prefix}_heatmap.png")
    plt.show()




# ====================== MAIN ======================

# === USER SETTINGS ===
FILE_PATH = '/Users/Hoan.Nguyen/ComBio/DebLab/ipi_ab_gpc.xlsx'
SEQUENCE_COLUMN = 'CDR3'                    # 'CDR3' or 'HSEQ'
DISTANCE_METHOD = 'levenshtein'             # 'levenshtein' or 'blosum62'
GROUP_BY_ANTIGEN = False

# === NEW: FLEXIBLE LABELS ===
# Write column names separated by colons to combine multiple columns into one label.
# Examples:
# 'BARCODE'
# 'BARCODE:SPR_Annot'
# 'SPR_Annot:BARCODE'
# 'antigen:BARCODE:SPR_Annot'

LABEL_COLUMNS= 'BARCODE'
#LABEL_COLUMNS = 'BARCODE:CDR3'   #← useful to see actual sequence
#LABEL_COLUMNS = 'BARCODE:SPR_Annot'         #← useful to see binding annotation

LABEL_SEPARATOR = ':'                       # Change to '_' or '-' if you prefer

COLOR_THRESHOLD = 2 if DISTANCE_METHOD == 'levenshtein' else 0.3

# === Load data ===
df = pd.read_excel(FILE_PATH)

## ALEX!!!
# you can uncomment or modify this line to limit number of sequences for testing
df =df[1:10]


# ALEX!!!
# you can change or remove or add more required columns here
required_cols = ['BARCODE', 'CDR3', 'HSEQ', 'SPR_Annot', 'antigen']
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Parse label columns
label_columns_list = [col.strip() for col in LABEL_COLUMNS.split(':') if col.strip()]

# Validate columns exist
missing_labels = [c for c in label_columns_list if c not in df.columns]
if missing_labels:
    raise ValueError(f"Label columns not found in file: {missing_labels}")

print(f"Loaded {len(df)} sequences")
print(f"Label format: {' : '.join(label_columns_list)}")

# === Run analysis ===
if GROUP_BY_ANTIGEN:
    for ag in df['antigen'].dropna().unique():
        subset = df[df['antigen'] == ag].reset_index(drop=True)
        if len(subset) < 2:
            print(f"Skipping antigen '{ag}' (only {len(subset)} sequence)")
            continue
        title_prefix = f"Antigen: {ag}"
        safe_ag = str(ag).replace(' ', '_').replace('/', '-').replace('\\', '-')
        save_prefix = f"cluster_{SEQUENCE_COLUMN}_{DISTANCE_METHOD}_{safe_ag}"
        
        plot_dendrogram_and_heatmap(
            subset, sequence_col=SEQUENCE_COLUMN,
            title_prefix=title_prefix,
            distance_method=DISTANCE_METHOD,
            color_threshold=COLOR_THRESHOLD,
            label_columns=label_columns_list,
            label_separator=LABEL_SEPARATOR,
            save_prefix=save_prefix
        )
else:
    title_prefix = "All Sequence Data"
    save_prefix = f"cluster_{SEQUENCE_COLUMN}_{DISTANCE_METHOD}_all"
    plot_dendrogram_and_heatmap(
        df, SEQUENCE_COLUMN, title_prefix=title_prefix,
        distance_method=DISTANCE_METHOD, color_threshold=COLOR_THRESHOLD,
        label_columns=label_columns_list, label_separator=LABEL_SEPARATOR,
        save_prefix=save_prefix
    )

print("\nYour beautiful tree Analysis complete!")


