# IBD regulatory networks

M2 BIM project (RESYS, Biological Networks course, 2026–2027)

In this project we compare gene regulatory networks between inflamed and non-inflamed intestinal tissue from patients with IBD (Crohn's disease and ulcerative colitis). We focus on genes from the TLR, NF-κB and JAK-STAT pathways, and we build the networks with PC and MIIC.
Authors: [Mouna Jegham](https://github.com/Mounaj3) and [Chloé Dubray](https://github.com/chloychlo)

## Data
We use the MSCCR cohort from GEO (GSE193677): bulk RNA-seq of about 2,500 intestinal biopsies.

Before starting, create a data/ folder and put these files in it:

The count matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193677/suppl/GSE193677_MSCCR_Biopsy_counts.txt.gz
The series matrix (metadata) and the gene annotation file (Human.GRCh38.p13.annot.tsv), both available here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193677
You don't need to unzip the count matrix, R can read it directly.
## What the code does

Run the scripts in this order:

1. `load_data.R`: clean the metadata, keep IBD patients only, and select the samples:
   all 712 inflamed biopsies, and 712 non-inflamed biopsies with the same disease × region
   distribution (stratified resampling, seed 123)
2. `data_preprocessing.R`: keep protein-coding genes, remove mitochondrial, ribosomal,
   constant and low-expressed genes, normalise (TMM + log2-CPM), and keep the 3,000 most
   variable genes 
3. `EDA.R`: explore the data: variance distribution, sample correlation and expression
   heatmaps, PCA (BioNERO)
5. `feature_selection__prep_data.R`: compute the mutual information between the 3,000 genes
   (infotheo, 10 equal-frequency bins), rank genes by mean MI and keep the top 40. Add 12
   regulators from the literature (Johnson et al. 2024, Argmann et al. 2023), normalised from
   the raw counts. Build the final MIIC input: 52 genes + 6 clinical variables
6. Build  networks with the PC and MIIC algorithms
7.Explore and visualise the results, compare the two networks, and run a pathway enrichment

The MI computation takes about 4.5 hours on 3,000 genes. The matrix is saved in
`results/MI_matrix.csv`, so you only need to compute it once.

## Key design choices

- Region is the main source of variance in this cohort. It is handled by stratified
  sampling and kept as a node in the network
- Clinical variables are network nodes: Disease, Region, Sex, endoscopic severity,
  CRP and histology score. MIIC handles categorical variables and missing values.
- Normalisation is done after the count-based filters

## Requirements

```r
install.packages(c("data.table", "dplyr", "ggplot2", "matrixStats",
                   "infotheo", "miic", "igraph", "circlize"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "edgeR", "BioNERO", "ComplexHeatmap"))
```

## References

- Argmann C. et al. Biopsy and blood-based molecular biomarker of inflammation in IBD.
  *Gut* 2023;72:1271–1287. (dataset)
- Johnson J.L. et al. Gene coexpression networks reveal a broad role for lncRNAs in
  inflammatory bowel disease. *JCI Insight* 2024;9:e168988.

