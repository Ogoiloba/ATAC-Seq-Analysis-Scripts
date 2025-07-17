#!/bin/bash

# Author: Ogo Iloba
# Description: Annotates ATAC-seq peaks with genomic features in Leucaena trichandra
# Requirements: BEDTools, awk, grep, sed
# Input: GFF3 annotation file and MACS2 peak summits
# Output: BED files of peaks in genomic regions and summary counts

#------------------------------#
# 1. Prepare Annotation Input  #
#------------------------------#

# Extract relevant feature types
grep -v "#" Ltri1.0_gene_models.gff3 | cut -f3 | sort | uniq

# Extract key fields: contig, feature type, start, end
cut -f1,3,4,5 Ltri1.0_gene_models.gff3 > annotation_raw.bed

# Keep only Contigs 0–27
awk '$1 ~ /^Contig(2[0-7]|[0-1]?[0-9])$/' annotation_raw.bed > annotation_filtered.bed

# Format: contig, start, end, feature
awk '{print $1"\t"$3"\t"$4"\t"$2}' annotation_filtered.bed > annotation_formatted.bed

#-----------------------------#
# 2. Generate Genomic Regions #
#-----------------------------#

# Exons
awk '$4 == "exon"' annotation_formatted.bed > exons.bed

# Gene coordinates with strand info
awk '$3 == "gene"' Ltri1.0_gene_models.gff3 | \
awk '{print $1 "\t" $4-1 "\t" $5 "\t" $7}' > genes_with_strand.bed

# Keep Contigs 0–27
awk '$1 ~ /^Contig(2[0-7]|[0-1]?[0-9])$/' genes_with_strand.bed > genes_filtered.bed

# Upstream (1 kb)
awk '{if ($4 == "+") print $1, ($2 - 1000 > 0 ? $2 - 1000 : 0), $2, $4; \
      else print $1, $3, $3 + 1000, $4}' OFS='\t' genes_filtered.bed > upstream.bed

# Downstream (1 kb)
awk '{if ($4 == "+") print $1, $3, $3 + 1000, $4; \
      else print $1, ($2 - 1000 > 0 ? $2 - 1000 : 0), $2, $4}' OFS='\t' genes_filtered.bed > downstream.bed

#------------------------#
# 3. Generate Introns    #
#------------------------#

# Filter full GFF3 for Contigs 0–27
awk '$1 ~ /^Contig([0-9]|1[0-9]|2[0-7])$/' Ltri1.0_gene_models.gff3 > annotation_contigs.gff3

# Extract exons and UTRs
awk '$3 == "exon" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" \
     { print $1 "\t" $4-1 "\t" $5 "\t" $9 }' annotation_contigs.gff3 > exons_utrs.bed

# Merge overlapping regions
sort -k1,1 -k2,2n exons_utrs.bed | bedtools merge -i stdin -c 4 -o collapse > merged_exons_utrs.bed

# Extract gene regions
awk '$3 == "gene" { print $1 "\t" $4-1 "\t" $5 "\t" $9 }' annotation_contigs.gff3 > genes.bed

# Subtract to get introns
bedtools subtract -a genes.bed -b merged_exons_utrs.bed > raw_introns.bed
awk '{print $1, $2, $3}' raw_introns.bed | sed 's/ \+/\t/g' > introns.bed

#-------------------------------------#
# 4. Prepare Peaks and Intersections  #
#-------------------------------------#

# From peak calling folder (e.g., C+E dataset)
cut -f1,2,3 opt_macs2_summits.bed > peaks.bed

# Filter peaks to Contigs 0–27
awk '$1 ~ /^Contig(2[0-7]|[0-1]?[0-9])$/' peaks.bed > filtered_peaks.bed

# Intersect peaks with each genomic feature
bedtools intersect -a filtered_peaks.bed -b exons.bed -wa > peaks_exons.bed
bedtools intersect -a filtered_peaks.bed -b introns.bed -wa > peaks_introns.bed
bedtools intersect -a filtered_peaks.bed -b upstream.bed -wa > peaks_upstream.bed
bedtools intersect -a filtered_peaks.bed -b downstream.bed -wa > peaks_downstream.bed

#-------------------------------------#
# 5. Calculate Counts & Intergenic    #
#-------------------------------------#

# Count peaks
total=$(wc -l < filtered_peaks.bed)
exon=$(wc -l < peaks_exons.bed)
intron=$(wc -l < peaks_introns.bed)
upstream=$(wc -l < peaks_upstream.bed)
downstream=$(wc -l < peaks_downstream.bed)

# Intergenic = total - sum of other categories
genic_total=$((exon + intron + upstream + downstream))
intergenic=$((total - genic_total))

echo "===== Peak Distribution Summary ====="
echo "Total Peaks: $total"
echo "Exonic Peaks: $exon"
echo "Intronic Peaks: $intron"
echo "Upstream Peaks: $upstream"
echo "Downstream Peaks: $downstream"
echo "Intergenic Peaks (inferred): $intergenic"

#-------------------------------------#
# 6. (Optional) Per-Chromosome Report #
#-------------------------------------#
# Run additional script (if available)
# ./Peakanno_per_chr_CE.sh
# Outputs: peak_annotation_table.txt, peak_annotation_percentages.txt

