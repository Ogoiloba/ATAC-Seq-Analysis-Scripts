#!/usr/bin/env python

"""
Author: Ogo Iloba
Script: Filter_genes_by_rpk.py
Description: From TPM-filtered gene pairs, retain only those where at least one gene has RPK >= 20.
Input: Chr0_21_filtered_by_tpm.tsv
Output: Final_Filtered_Data_021.tsv
"""

import pandas as pd

# Load the filtered TPM data
input_file = "Chr0_21_filtered_by_tpm.tsv"
df = pd.read_csv(input_file, delimiter='\t')

# Convert RPK columns to numeric
for col in ['RPK1', 'RPK2']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Drop rows with missing RPK values
df.dropna(subset=['RPK1', 'RPK2'], inplace=True)

# Keep rows where either RPK1 or RPK2 is >= 20
df = df[(df['RPK1'] >= 20) | (df['RPK2'] >= 20)]

# Save the final filtered data
df.to_csv("Final_Filtered_Data_021.tsv", sep='\t', index=False)

# Debug output
print("Final filtered data (first 10 rows):")
print(df.head(10))
print(f"\nNumber of final filtered pairs: {df.shape[0]}")
