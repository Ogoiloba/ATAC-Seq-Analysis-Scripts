### ATAC-Seq: Trimming and Mapping Pipeline with FRiP Score – Dataset A

This script performs a complete ATAC-seq processing pipeline using Dataset A, from raw FASTQ files to peak calling and FRiP (Fraction of Reads in Peaks) score calculation. This pipeline was used for the other libraries C,D & E. 

## Steps:

1. Trimming sequencing adapters and low-quality reads using Trim Galore
2. Alignment to the Leucaena genome using Burrow-Wheeler Aligner (BWA)
3. Filtering low-quality reads and removing PCR duplicates with SAMtools
4. ATAC-seq Read Shifting using deepTools `alignmentSieve`
5. Peak Calling with MACS2 and merging nearby peaks (<10 bp apart) using BEDTools
6. FRiP Calculation- proportion of ATAC-seq reads found in peaks

## Input files
- Paired-end FASTQ files:
/media/CLG2018/ogoiloba/ATAC-seq/Xuan_2_NS_S30_L002_R1_001.fastq
/media/CLG2018/ogoiloba/ATAC-seq/Xuan_2_NS_S30_L002_R2_001.fastq

 Reference genome:
/media/CLG2018/ogoiloba/ATAC-seq/Ltri0.3_genome.fasta

## Outputs
- Trimmed and aligned BAM files
- MACS2 peak files (`.xls`, `.bed`)
- Merged peak file (`DatasetA_peaks_Lengths.bed`)
- FRiP score report:
- Total mapped reads
- Reads in peaks
- Fraction and percentage in peaks
