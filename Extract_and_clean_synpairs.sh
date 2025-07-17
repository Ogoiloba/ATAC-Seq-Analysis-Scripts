#!/bin/bash
# Author: Ogo Iloba
# Project: ATAC-seq - Peak Length and Gene Expression Analysis (Datasets C & E)
# Description: Step 1 - Extract syntenic gene pairs from SynMap data, remove duplicates, and map them to original synteny blocks.

#---------------------------#
# Step 1: Extract syntenic gene pairs with block info
#---------------------------#
awk '{if ($1 ~ /^#/) print $1, $3, $4, $5, $6; else print $5, $17}' Tab_Genocord_grpanalysis2.txt > Synteny_with_GenesOnly.txt

# Keep only Contig0 - Contig21 pairs
awk '/^#/ {if ($2 ~ /_Contig0$/ && $3 ~ /_Contig21$/) p=1; else p=0} p' Synteny_with_GenesOnly.txt > 021_Synteny_with_GenesOnly.txt

# Remove headers and convert to tab-delimited
grep -v "#" 021_Synteny_with_GenesOnly.txt | sed 's/ \+/\t/g' > noheader_021_Synteny_with_GenesOnly.txt

# Extract headers separately
grep "^#" 021_Synteny_with_GenesOnly.txt > header_names_021.txt

# Sort and remove duplicated rows
sort noheader_021_Synteny_with_GenesOnly.txt | uniq > sorted_noheader_021_Synteny_with_GenesOnly.txt

#---------------------------#
# Step 2: Remove duplicate gene pairs (many-to-one)
#---------------------------#
awk '{count1[$1]++; count2[$2]++; pairs[$1,$2]++} \
     END {for (key in pairs) {split(key, arr, SUBSEP); \
     if ((count1[arr[1]] == 1 && count2[arr[2]] == 1) || pairs[key] > 1) \
     print arr[1], arr[2]}}' sorted_noheader_021_Synteny_with_GenesOnly.txt > Duplicates_removed_021.txt

# Convert to tab-delimited
sed 's/ \+/\t/g' Duplicates_removed_021.txt > duplicate_genes_removed_021.txt

#---------------------------#
# Step 3: Map gene pairs back to synteny block headers
#---------------------------#
input_file="duplicate_genes_removed_021.txt"
original_data="merged_GroupAB0_21"
output_file="mapped_duplicate_genes_removed_021.txt"
> "$output_file"

while IFS= read -r line; do
    if [[ "$line" =~ ^# ]]; then
        current_block="$line"
        echo "$line" >> "$output_file"
    else
        gene1=$(echo "$line" | awk '{print $1}')
        gene2=$(echo "$line" | awk '{print $2}')
        if grep -q -E "^$gene1\s+$gene2$" "$input_file"; then
            echo "$line" >> "$output_file"
            sed -i "/^$gene1\s\+$gene2$/d" "$input_file"
        fi
    fi
done < "$original_data"

# Extract final gene pair list with synteny headers cleaned
awk '{if ($1 ~ /^#/) print $1, $2, $3, $4; else print $4, $13}' mapped_duplicate_genes_removed_021.txt | sed 's/\.1//g' > Mapped021_nodups_Synteny_with_GenesOnly.txt
