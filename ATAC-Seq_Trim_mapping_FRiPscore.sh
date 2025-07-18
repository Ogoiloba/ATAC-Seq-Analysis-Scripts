#!/bin/bash
# Author: Ogo Iloba
# Description: Full ATAC-seq pipeline for library A with actual file paths. This same script was used for the other datasets.
# Input: Raw paired-end FASTQ files
# Output: BAM files, peaks, and FRiP score


#SETUP          

# Raw data
RAW_R1="/media/CLG2018/ogoiloba/ATAC-seq/Xuan_2_NS_S30_L002_R1_001.fastq"
RAW_R2="/media/CLG2018/ogoiloba/ATAC-seq/Xuan_2_NS_S30_L002_R2_001.fastq"

# Genome reference
GENOME="/media/CLG2018/ogoiloba/ATAC-seq/Ltri0.3_genome.fasta"

# Output directories
TRIM_DIR="/media/CLG2018/ogoiloba/ATAC-seq/01_trim_fastq"
WORK_DIR="/media/CLG2018/ogoiloba/ATAC-seq"

mkdir -p "$TRIM_DIR"


#1. TRIMMING        

trim_galore --paired --fastqc --gzip --length 15 -q 20 -j 1 \
  --basename trim -o "$TRIM_DIR" "$RAW_R1" "$RAW_R2"

# Move trimmed files back to working dir
cp "$TRIM_DIR/trim_val_1.fq.gz" "$WORK_DIR"
cp "$TRIM_DIR/trim_val_2.fq.gz" "$WORK_DIR"

cd "$WORK_DIR"
gzip -d trim_val_1.fq.gz
gzip -d trim_val_2.fq.gz


# 2. ALIGNMENT      
bwa mem "$GENOME" trim_val_1.fq trim_val_2.fq -t 24 > datasetA_BWA_mapped_reads.sam

samtools view -bS datasetA_BWA_mapped_reads.sam > datasetA_BWA_mapped_reads.bam

# Alignment stats
samtools flagstat datasetA_BWA_mapped_reads.bam > datasetA_flagstat.txt


# 3. REMOVE PCR DUPLICATES    
samtools view -h -q 30 datasetA_BWA_mapped_reads.bam | \
  samtools view -b -o datasetA_BWA_mapped_qreads.bam

samtools sort -n datasetA_BWA_mapped_qreads.bam -o datasetA_BWA_qreads_sorted.bam
samtools fixmate -m datasetA_BWA_qreads_sorted.bam datasetA_BWA_qreads_fixmate.bam
samtools sort -o datasetA_BWA_qreads_positionsort.bam datasetA_BWA_qreads_fixmate.bam
samtools markdup -r -s datasetA_BWA_qreads_positionsort.bam datasetA_BWA_qreads_markdup.bam


# 4. ATAC SHIFT        
samtools index datasetA_BWA_qreads_markdup.bam
alignmentSieve --ATACshift --bam datasetA_BWA_qreads_markdup.bam -o datasetA_BWA_qreads_markdup_shift.bam

cp datasetA_BWA_qreads_markdup_shift.bam "$TRIM_DIR"


# 5. FINAL SORTING + STATS     
cd "$TRIM_DIR"
samtools sort -O bam -o datasetA_BWA_qreads_markdup_shiftSort.bam datasetA_BWA_qreads_markdup_shift.bam
samtools index datasetA_BWA_qreads_markdup_shiftSort.bam
samtools sort -n -O bam -o datasetA_BWA_qreads_markdup_shiftSortN.bam datasetA_BWA_qreads_markdup_shiftSort.bam

samtools stats datasetA_BWA_qreads_markdup.bam > datasetA_BWA_qreads_markdup_summary.txt
samtools stats datasetA_BWA_qreads_markdup_shiftSortN.bam > datasetA_BWA_qreads_markdup_shiftSortN_summary.txt


# 6. PEAK CALLING & MERGE PEAKS < 10bp APART    
mkdir -p peak_calling
cd peak_calling

macs2 callpeak --cutoff-analysis \
  -t ../datasetA_BWA_qreads_markdup_shiftSortN.bam \
  -f BAM -g 8.5e8 -q 0.05 --bdg -n opt_macs2

grep -v "#" opt_macs2_peaks.xls | cut -f1,2,3,4,5,10 > datasetA_macs2_peaks.xls
grep -v "chr" datasetA_macs2_peaks.xls > DatasetA_Macs2_Peaks.xls

bedtools merge -i DatasetA_Macs2_Peaks.xls -d 10 -c 4,5,6 -o sum,collapse,collapse > DatasetA_peaks_Lengths.bed


# 7. FRiP SCORE CALCULATION     
bedtools bamtobed -i ../datasetA_BWA_qreads_markdup_shiftSortN.bam > A_Atac_reads.bed
sed 's/ \+/\t/g' A_Atac_reads.bed > DatasetA_Atac_reads.bed

bedtools intersect -a DatasetA_Atac_reads.bed -b DatasetA_peaks_Lengths.bed -u > DatasetA_Atacreads_in_peaks.txt

# Calculate FRiP
total_reads=$(wc -l < DatasetA_Atac_reads.bed)
reads_in_peaks=$(wc -l < DatasetA_Atacreads_in_peaks.txt)
fraction=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)
percentage=$(echo "$fraction * 100" | bc)

echo "Dataset A FRiP Report"
echo "---------------------"
echo "Total mapped reads: $total_reads"
echo "Reads found in peaks: $reads_in_peaks"
echo "Fraction of reads in peaks: $fraction"
echo "Percentage of reads in peaks: $percentage%"


