# Bioinformatics Scripts

This repository contains Rmds, docs bash scipts create for BIOL 668. Practice with genome preprocessing, RNA-seq analysis, and microbiome analysis using QIIME2. Tools used include `fastp`, `Kaiju`, `QIIME2`, `DESeq2`, ...

---

# Folder Descriptions

# AlgsAndGenomes/
- Quality filtering for a metagenomic dataset.
- Used `fastp` to clean and trim FASTQ 
- Generated quality reports in HTML/JSON and assessed GC content and k-mer frequencies.
- Attempted to classify sequences 
- Explored data stability based on GC content and quality metrics.

# QIIME2/
- Microbiome diversity and taxonomy workflow using QIIME2.
- Imported and demultiplexed sequencing data 
- Denoised sequences with DADA2 
- Performed alpha and beta diversity analysis (Faith PD, Bray Curtis, UniFrac).
- Conducted differential abundance analysis of gut microbiome using ANCOM-BC.
- Created taxa barplots and Emperor plots to visualize diversity and sample clustering.
- Step-by-step commands documented in `QIIME2.txt` and results summarized in `QIIME2 part2.docx`.

# RNA_SEQ_LAB/
- Rmd files for RNA-seq differential expression analysis.
- Read alignment and count summarization 
- Used DESeq2 in R to compare gene expression between mutant and wild type samples.
- Included exploratory data analysis and visualizations.
---

