#!/bin/bash

# Author: Ogo Iloba
# Script: prepare_peak_overlap.sh
# Description: Extract gene coordinates for duplicates and identify and merge overlapping peaks overlapping with genes from datasets C and E.

#---------------------------#
# Step 1: Prepare gene lists from filtered file
#---------------------------#
cut -f1 Final_Filtered_Data_021.tsv | grep -v "GENE1" > homeolog_goi_0.txt
cut -f6 Final_Filtered_Data_021.tsv | grep -v "GENE2" > homeolog_goi_21.txt

#---------------------------#
# Step 2: Extract gene coordinates from annotated gene BED
#---------------------------#
awk 'BEGIN{FS=OFS="\t"} NR==FNR{genes[$1]; next} $2=="gene" && match($6, /ID=([^;]+)/, m) && (m[1] in genes) {
  print $1, $3, $4, m[1]
}' homeolog_goi_0.txt tabLtri_up_downstream_genes.bed > gene_cords_Contig0.bed

awk 'BEGIN{FS=OFS="\t"} NR==FNR{genes[$1]; next} $2=="gene" && match($6, /ID=([^;]+)/, m) && (m[1] in genes) {
  print $1, $3, $4, m[1]
}' homeolog_goi_21.txt tabLtri_up_downstream_genes.bed > gene_cords_Contig21.bed

#---------------------------#
# Step 3: Prepare peak files
#---------------------------#
grep -v "#" opt_macs2_peaks.xls | cut -f1,2,3,4,5,10 > macs2_peaks.xls
grep -v "chr" macs2_peaks.xls > Macs2_Peaks.xls

# Merge nearby peaks (<= 10 bp)
bedtools merge -i Macs2_Peaks.xls -d 10 -c 4,5,6 -o sum,collapse,collapse > Peaks_Lengths_CE.bed

#---------------------------#
# Step 4: Intersect peaks with gene coordinates
#---------------------------#
bedtools intersect -a gene_cords_Contig0.bed -b Peaks_Lengths_CE.bed -wa -wb > Contig0_overlaps_CE.bed
bedtools intersect -a gene_cords_Contig21.bed -b Peaks_Lengths_CE.bed -wa -wb > Contig21_overlaps_CE.bed

#---------------------------#
# Step 5: Merge peak info per gene
#---------------------------#
sort -k4,4 -k6,6n Contig0_overlaps_CE.bed > Sorted_Contig0_overlaps_CE.bed
sort -k4,4 -k6,6n Contig21_overlaps_CE.bed > Sorted_Contig21_overlaps_CE.bed

# Merge peaks per gene block (Contig0)
awk 'BEGIN { FS = OFS = "\t" }
{
  if ($4 != current_gene) {
    if (NR > 1) {
      print current_contig, gene_start, gene_stop, current_gene, peak_contig, peak_start, peak_stop, peak_sum, peak_summits, peak_ids;
    }
    current_gene = $4;
    current_contig = $1;
    gene_start = $2;
    gene_stop = $3;
    peak_contig = $5;
    peak_start = $6;
    peak_stop = $7;
    peak_sum = $8;
    peak_summits = $9;
    peak_ids = $10;
  } else {
    peak_stop = $7;
    peak_sum += $8;
    peak_summits = peak_summits "," $9;
    peak_ids = peak_ids "," $10;
  }
}
END {
  print current_contig, gene_start, gene_stop, current_gene, peak_contig, peak_start, peak_stop, peak_sum, peak_summits, peak_ids;
}' Sorted_Contig0_overlaps_CE.bed > final_Contig0_overlaps_CE.bed

# Repeat for Contig21 (can be refactored if needed)
awk 'BEGIN { FS = OFS = "\t" }
{
  if ($4 != current_gene) {
    if (NR > 1) {
      print current_contig, gene_start, gene_stop, current_gene, peak_contig, peak_start, peak_stop, peak_sum, peak_summits, peak_ids;
    }
    current_gene = $4;
    current_contig = $1;
    gene_start = $2;
    gene_stop = $3;
    peak_contig = $5;
    peak_start = $6;
    peak_stop = $7;
    peak_sum = $8;
    peak_summits = $9;
    peak_ids = $10;
  } else {
    peak_stop = $7;
    peak_sum += $8;
    peak_summits = peak_summits "," $9;
    peak_ids = peak_ids "," $10;
  }
}
END {
  print current_contig, gene_start, gene_stop, current_gene, peak_contig, peak_start, peak_stop, peak_sum, peak_summits, peak_ids;
}' Sorted_Contig21_overlaps_CE.bed > final_Contig21_overlaps_CE.bed
