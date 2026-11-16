# IBD regulatory networks

M2 BIM project (RESYS, Biological Networks course, 2026–2027)

In this project we compare gene regulatory networks between inflamed and non-inflamed
intestinal tissue from patients with IBD (Crohn's disease and ulcerative colitis).
We focus on genes from the TLR, NF-κB and JAK-STAT pathways, and we build the networks
with PC and MIIC.

Authors: [Mouna Jegham](github.com/Mounaj3) and [Chloé Dubray](github.com/chloychlo)

## Data

We use the MSCCR cohort from GEO (GSE193677): bulk RNA-seq of about 2,500 intestinal biopsies.

**Before starting**, create a `data/` folder and put these files in it:

- The count matrix:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193677/suppl/GSE193677_MSCCR_Biopsy_counts.txt.gz
- The series matrix (metadata) and the gene annotation file (`Human.GRCh38.p13.annot.tsv`),
  both available here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193677

You don't need to unzip the count matrix, R can read it directly.

## What the code does

1. Clean the metadata and select the samples: inflamed biopsies, and the same number of
   non-inflamed biopsies with the same disease and region distribution (stratified resampling).
2. Preprocess the counts: keep protein-coding genes, remove mitochondrial and ribosomal
   genes and low-expressed genes, then normalise (TMM + log-CPM).
3. Explore the data: PCA and correlation heatmaps (BioNERO).
4. Gene selection: genes of the TLR, NF-κB and JAK-STAT pathways, keeping all transcription
   factors and the most variable genes (MAD on region-corrected data).
5. Build two networks, one for inflamed tissue and one for non-inflamed tissue, with the
   PC and MIIC algorithms.
6. Explore and visualise the results, compare the two networks, and run a pathway enrichment.


## Description of the R scripts
* **01_load_data**: Load the GEO dataset GSE193677, extract and clean sample metadata, and select a balanced set of inflamed and non-inflamed IBD biopsy samples.
  * Inputs: _data/raw/GSE193677_series_matrix.txt_	(GEO series matrix containing sample metadata), _data/raw/GSE193677_MSCCR_Biopsy_counts.txt_	(Raw RNA-seq gene expression count matrix)
  * Outputs: _data/processed/counts_selected_samples.csv_	(Raw gene expression counts for selected samples), _data/processed/meta_selected_samples.csv_	(Cleaned metadata for selected samples)
* **02_data_preprocessing**: Filter genes based on their biotype, remove mitochondrial and ribosomal genes, normalise RNA-seq counts using TMM and log2-CPM, and retain the 3,000 most variable genes.
  * Inputs: _data/processed/counts_selected_samples.csv_ (Raw gene expression counts), _data/raw/Human.GRCh38.p13.annot.tsv_ (Gene annotation file)
  * Outputs: _data/processed/expr_clean_logcpm.csv_ (Filtered and normalised gene expression matrix), _data/processed/selected_genes_annotation.csv_ (Annotations of retained genes), _data/processed/preprocessed_data.rds_ (Combined preprocessed data and annotations)
* **03_network_inference_PC**: Select genes from inflammation-related Reactome pathways (including transcription factors and the 100 most variable pathway genes), infer separate gene networks using the PC algorithm for inflamed and non-inflamed samples, and export network objects, edge lists and visualisations.
  * Inputs: _data/processed/preprocessed_data.rds_ (Preprocessed expression matrix, gene annotations and sample metadata).
  * Outputs: _data/processed/selected_genes_annotation_PC.csv_ (Annotations of selected genes), _data/processed/selected_expression_matrices.rds_ (Expression matrices for both conditions), _results/networks/PC_network_ (Inferred PC network objects and edge lists), _results/figures/PC_network.png_ (Network visualisations).
* **04_network_analysis_PC**: Compare the inflamed and non-inflamed PC networks by analysing their topology, connectivity, shared and condition-specific interactions, and performing exploratory Reactome pathway enrichment on condition-specific genes.
  * Inputs: _results/networks/PC_network_inflamed_edges.csv_ and _results/networks/PC_network_noninflamed_edges.csv_ (Inferred network edge lists), _data/processed/selected_expression_matrices.rds_ (Selected gene annotations).
  * Outputs: _results/analysis/tables/_ (CSV files containing network statistics, node degrees, edge comparisons and Reactome enrichment results), _results/analysis/networks/_ (RDS network and enrichment objects), _results/analysis/figures/_ (Network visualisations and Reactome enrichment plots).

## Key design choices

- Normalisation is done after the count-based filters, on both groups together, so the two
  networks share the same scale.
- Region is the main source of variance in this cohort. It is handled
  by stratified sampling.

## Requirements

```r
# CRAN packages
install.packages(c(
  "data.table",
  "dplyr",
  "matrixStats",
  "msigdbr",
  "igraph",
  "ggraph",
  "ggplot2",
  "ggrepel"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GEOquery",
  "BioNERO",
  "edgeR",
  "pcalg",
  "dorothea",
  "clusterProfiler"
))
```

## References

- Argmann C. et al. Biopsy and blood-based molecular biomarker of inflammation in IBD.
  *Gut* 2023;72:1271–1287. (dataset)
- Johnson J.L. et al. Gene coexpression networks reveal a broad role for lncRNAs in
  inflammatory bowel disease. *JCI Insight* 2024;9:e168988.
