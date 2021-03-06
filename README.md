# thesis-RNAseq
Analysis of B. saida full transcriptome data from raw reads to gene ontology

Steps outlined in this repository:
1. Demultiplexing
2. QA/QC
3. Adapter trimming with Bbduk
4. Genome indexing and alignment with STAR
5. Feature count with htseq-count
6. Verification of sample inter-comparability
7. Differential gene expression (DGE) analysis
   * identification of up and downregulated genes with DESeq2
   * conversion of genes to standardized UIDs with NCBI's eFetch Utility
   * annotation of biological processes and assignment of gene ontologies using FishEnrichr
   * visualization of GO pathways using simplifyEnrichment 

Steps 1-5 can be found in the [pipeline](pipeline.md), while inter-comparability verification is carried out in R (version 4.0.3 (2020-10-10)) through the approach described in [fish stats](fish_stats_markdown.md) and DGE is documented in the [DESeq2 markdown](DESeq2-markdown.md). For additional reading, including a deeper discussion of these analyses, check out my MSc. [thesis](https://hdl.handle.net/10037/21751).
