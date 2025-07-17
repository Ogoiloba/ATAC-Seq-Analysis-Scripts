#!/usr/bin/env python

"""
Author: Ogo Iloba
Script: Filter_genes_by_tpm.py
Description: Filters homeologous gene pairs to retain those with at least a 2-fold difference in TPM.
Input: Merged_TPM_021_file.txt (tab-delimited with headers)
Output: Chr0_21_filtered_by_tpm.tsv
"""

import pandas as pd
import numpy as np

# Load the dataset
input_file = "Merged_TPM_021_file.txt"
df = pd.read_csv(input_file, delimiter='\t', comment='#')

# Ensure TPM columns are numeric
for col in ['TPM1', 'TPM2']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Drop rows with missing or zero TPM values
df.dropna(subset=['TPM1', 'TPM2'], inplace=True)
df = df[(df['TPM1'] != 0) & (df['TPM2'] != 0)]

# Filter for 2-fold difference in TPM
mask = (df['TPM1'] >= 2 * df['TPM2']) | (df['TPM2'] >= 2 * df['TPM1'])
filtered = df[mask]

# Save result
filtered.to_csv("Chr0_21_filtered_by_tpm.tsv", sep='\t', index=False)

# Debug output
print("Filtered data (first 10 rows):")
print(filtered.head(10))
print(f"\nNumber of filtered pairs: {filtered.shape[0]}")
