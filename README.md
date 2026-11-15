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

## Key design choices

- Normalisation is done after the count-based filters, on both groups together, so the two
  networks share the same scale.
- Region is the main source of variance in this cohort. It is handled
  by stratified sampling.

## Requirements

```r
install.packages(c("data.table", "dplyr", "tidyr", "ggplot2", "matrixStats",
                   "reshape2", "miic", "igraph"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "edgeR", "limma", "BioNERO",
                       "KEGGREST", "decoupleR", "pcalg"))
```

## References

- Argmann C. et al. Biopsy and blood-based molecular biomarker of inflammation in IBD.
  *Gut* 2023;72:1271–1287. (dataset)
- Johnson J.L. et al. Gene coexpression networks reveal a broad role for lncRNAs in
  inflammatory bowel disease. *JCI Insight* 2024;9:e168988.
