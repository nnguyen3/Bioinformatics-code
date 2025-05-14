# Bioinformatics Scripts

This repository contains Rmds, docs bash scipts create for BIOL 668 at SDSU. Practice with genome preprocessing, RNA-seq analysis, and microbiome analysis using QIIME2. Tools used include `fastp`, `Kaiju`, `QIIME2`, `DESeq2`, ...

---

# Folder Descriptions

### AlgsAndGenomes/
 # Required file for AlgsAndGenomes/
- out.insub732_2_R1_fastp.fastq (forward R1, paired-end read)
- insub732_2_R2.fastp.fastq.gz (reverse R2, paired-end read)
- reads_1.fq (test file)
- reads_2.fq (test file) 
- viruses database provided by the tutorial
- To create conda environment and activate fastp
```
conda creaete -n fastp
conda  activate fastp
conda install -c bioconda fastp

```
-To insatll Kaiju
```
conda install -c bioconda kaiju
```
- Quality filtering for a metagenomic dataset.
- Used `fastp` to clean and trim FASTQ 
- Generated quality reports in HTML/JSON and assessed GC content and k-mer frequencies.
- Attempted to classify sequences 
- Explored data stability based on GC content and quality metrics.

### QIIME2/
# Required file for QIIME2/
- https://data.qiime2.org/2024.10/tutorials/moving-pictures/sample_metadata.tsv
- https://data.qiime2.org/2024.10/tutorials/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz
- https://data.qiime2.org/2024.10/tutorials/moving-pictures/emp-single-end-sequences/sequences.fastq.gz
- Use tutorial to working with QIIME2 tool : https://docs.qiime2.org/2024.10/tutorials/moving-pictures/.
- To install QIIME2 : https://docs.qiime2.org/2024.10/install/.
- Microbiome diversity and taxonomy workflow using QIIME2.
- Imported and demultiplexed sequencing data 
- Denoised sequences with DADA2 
- Performed alpha and beta diversity analysis (Faith PD, Bray Curtis, UniFrac).
- Conducted differential abundance analysis of gut microbiome using ANCOM-BC.
- Created taxa barplots and Emperor plots to visualize diversity and sample clustering.
- Step-by-step commands documented in `QIIME2.txt` and results summarized in `QIIME2 part2.docx`.

### RNA_SEQ_LAB/
# Required file for RNA_SEQ_LAB/
- rna_counts_data.csv
- rna_map_update copy.csv
- pbmc3k_filtered_gene_bc_matrices.tar.gz
- part 1 is to follow the Bioconductor project's tutorial (https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html. )
- Use Glimma, limma, and edgeR in R to filter, organize, plot, and investigate the data
- part 2 is use tutorial and use DESeq2 to perform Differential Gene Expression (DGE) analysis, between a treatment and control group
- to install DESeq2
```
# install DESeq2
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

# install tidyverse (ggplot2 is in tidyverse)
install.packages("tidyverse")

```
- part 3 is to analyze single cell RNA with Seurat , tutorial (https://satijalab.org/seurat/.)
- to instal Seurat
```
install.packages('Seurat')
library(Seurat)

```
- Rmd files for RNA-seq differential expression analysis.
- Read alignment and count summarization 
- Used DESeq2 in R to compare gene expression between mutant and wild type samples.
- Included exploratory data analysis and visualizations.
---

