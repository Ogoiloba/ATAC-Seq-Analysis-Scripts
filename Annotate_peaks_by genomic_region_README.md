# ATAC-Seq Peak Annotation Script

This script annotates ATAC-seq peak locations (from MACS2) by intersecting them with genomic regions including exons, introns, upstream, downstream, and intergenic regions in Leucaena trichandra.

## Steps:

1. Parse and format GFF3 annotations
2. Extract genomic regions:
   - Exons
   - Introns (by subtracting exons from genes)
   - Upstream/downstream regions (1 kb)
3. Filter peaks to major contigs (0–27) (We did this because the gff3 file contains more than 27 contigs, but the actual chromomosomes are 0- 27
4. Intersect peaks with each feature
5. Count peaks per category and infer intergenic

## Inputs

 – genome annotation file - `Ltri1.0_gene_models.gff3`
-  ATAC-seq peaks from MACS2 -`opt_macs2_summits.bed` 

## Outputs

- Bed files of the peaks in the different genomic regions.
- A summary printed to screen:
  - Total peaks
  - Number in each category
  - Intergenic peaks (calculated)

