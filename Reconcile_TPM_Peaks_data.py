#!/usr/bin/env python

"""
Author: Ogo Iloba
Script: testreconcile_TPM_Peak_Values.py
Description: Merge TPM/RPK data and gene-level peak information into a single table for plotting.
Input: TPM_RPK_Homeologs_CE.txt, PeaksInContig0_genes_CE.txt, PeaksInContig21_genes_CE.txt, gene_cords_Contig0.bed, gene_cords_Contig21.bed
Output: reconciled_data_021_CE.txt
"""

import pandas as pd

# Load homeologs data
homeologs_df = pd.read_csv('TPM_RPK_Homeologs_CE.txt', sep='\t', header=None)
homeologs_df.columns = ['GENE1', 'RPK1', 'TPM1', 'GENE2', 'RPK2', 'TPM2']

# Load coordinates for Contig0 and Contig21 genes
coords_contig0 = pd.read_csv('gene_cords_Contig0.bed', sep='\t', header=None, names=['Contig', 'GeneStart', 'GeneStop', 'GeneID'])
coords_contig21 = pd.read_csv('gene_cords_Contig21.bed', sep='\t', header=None, names=['Contig', 'GeneStart', 'GeneStop', 'GeneID'])

# Load peak data per gene
peaks_contig0 = pd.read_csv('PeaksInContig0_genes_CE.txt', sep='\t', header=None)
peaks_contig0.columns = ['GeneStartPosition', 'GeneID', 'Peak_length', 'peak_ID']

peaks_contig21 = pd.read_csv('PeaksInContig21_genes_CE.txt', sep='\t', header=None)
peaks_contig21.columns = ['GeneStartPosition', 'GeneID', 'Peak_length', 'peak_ID']

# Function to merge homeologs with coordinates and peaks
def merge_with_data(df, coords, peaks, gene_col, suffix):
    merged = pd.merge(df, coords, left_on=gene_col, right_on='GeneID', how='left')
    merged.rename(columns={'GeneStart': f'GeneStart{suffix}'}, inplace=True)

    peak_summary = peaks.groupby('GeneID').agg({
        'Peak_length': lambda x: ','.join(map(str, x)),
        'peak_ID': lambda x: ','.join(map(str, x))
    }).reset_index()

    merged = pd.merge(merged, peak_summary, on='GeneID', how='left')
    merged.drop(columns=['GeneID', 'Contig', 'GeneStop'], inplace=True)

    merged.rename(columns={
        'Peak_length': f'Peak_length{suffix}',
        'peak_ID': f'peak_ID{suffix}'
    }, inplace=True)

    for col in [f'Peak_length{suffix}', f'peak_ID{suffix}']:
        merged[col].fillna('0', inplace=True)

    return merged

# Merge data for both GENE1 and GENE2
merged = merge_with_data(homeologs_df, coords_contig0, peaks_contig0, 'GENE1', '1')
merged = merge_with_data(merged, coords_contig21, peaks_contig21, 'GENE2', '2')

# Ensure consistent column order
columns_order = ['GENE1', 'RPK1', 'TPM1', 'GeneStart1', 'Peak_length1', 'peak_ID1',
                 'GENE2', 'RPK2', 'TPM2', 'GeneStart2', 'Peak_length2', 'peak_ID2']
for col in columns_order:
    if col not in merged.columns:
        merged[col] = '0'

# Save merged output
merged.to_csv('reconciled_data_021_CE.txt', sep='\t', index=False)

print("Reconciled data saved to reconciled_data_021_CE.txt")
