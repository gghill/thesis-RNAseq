# thesis-RNAseq
Analysis of B. saida full transcriptome data from raw reads to gene ontology

`Pipeline documentation in progress...[###  `   ` ` ` `  `  ]`

Steps outlined in this repository:
1. Demultiplexing
2. QA/QC
3. Adapter trimming with Bbduk
4. Genome indexing and alignment with STAR
5. Feature count with htseq-count
6. Differential gene expression analysis
   * identification of up and downregulated genes with DESeq2
   * conversion of genes to standardized UIDs with NCBI's eFetch Utility
   * annotation of biological processes and assignment of gene ontologies using FishEnrichr
   * visualization of GO pathways using simplifyEnrichment 
