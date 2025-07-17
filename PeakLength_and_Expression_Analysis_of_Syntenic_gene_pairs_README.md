# Peak Length and Expression Analysis of Syntenic Gene Pairs in L.trichandra.

This set of scripts performs an analysis on syntenic gene pairs across Contigs in Leucaena trichandra (Contigs 0 and 21 are used as examples here)  with their ATAC-seq peaks and expression values (TPM, RPK) The goal is to identify relationships between chromatin accessibility (measured by the length of peaks) and gene expression across duplicated gene blocks.The information from this was used to make graphs to visualize this relationship and also carryout statistical tests.

## Summary of Scripts
# 1. `Extract_and_clean_synpairs.sh`
- Extracts syntenic gene pairs between Contig0 and Contig21 from CoGe SynMap output.
- Removes many-to-one duplicate pairs to retain unique mappings.
- Re-maps clean gene pairs back to their original synteny block headers.
- Final output: `Mapped021_nodups_Synteny_with_GenesOnly.txt`

# 2. `Filter_by_TPM.py` and `Filter_by_RPK.py`
- Takes expression tables and filters the syntenic gene pairs using:
  - TPM (Transcripts Per Million) and
  - RPK (Reads Per Kilobase)thresholds
- Produces a cleaned file (`Final_Filtered_Data_021.tsv`) of expressed gene pairs for downstream peak analysis.

# 3. `Merge_Peaks.sh`
- Extracts peak data from ATAC-seq datasets.
- Merges overlapping peaks within 10bp and intersects them with gene coordinates.
- Collapses peak signals for each gene and creates final peak summary files:
  - `final_Contig0_overlaps_CE.bed`
  - `final_Contig21_overlaps_CE.bed`

# 4. `Reconcile_TPM_Peaks_data.py`
- Integrates the filtered expression data with merged peak data.
- Aligns gene IDs from Contig0 and Contig21 and combines peak and expression information into a single table.
- Final output: a unified dataset for comparing chromatin accessibility and gene expression between homeologs.

