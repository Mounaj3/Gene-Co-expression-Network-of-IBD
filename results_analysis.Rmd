# Results analysis

This document presents and guides the interpretation of the gene
association networks inferred from inflamed and non-inflamed IBD biopsy
samples using the PC algorithm.

The analysis is exploratory: the networks represent
conditional-dependence relationships inferred from the expression data.
The PC algorithm estimates a **CPDAG** (completed partially directed
acyclic graph).


------------------------------------------------------------------------

## 1. Network inference

Separate PC networks were inferred for inflamed and non-inflamed samples
using the same selected gene set. This shared gene set makes the two
inferred networks more directly comparable: differences in their edge
sets are not simply due to comparing different sets of candidate genes.

The selected genes were drawn from inflammation-related Reactome
pathways (including TLR, NF-κB and JAK-STAT-related pathways), with
transcription factors and highly variable pathway genes included in the
selection. Region-corrected variability was used for ranking genes,
whereas the original log-CPM expression values were used for network
inference.

### 1.1 Inflamed samples

![PC network inferred from inflamed
samples](results/figures/PC_network_inflamed.png)

Most connected nodes in the inflamed network:

| gene     | degree  | 
| -------- | ------- |
| BTK      | 8       |
| LCP1     | 7       |
| IL10     | 7       |
| LY96     | 7       |
| DUSP4    | 6       |
| NOD2     | 6       |


It appears from the inferred network that the BTK gene has the highest degree with 8 connections.

**Interpretation**
-   Are there visibly distinct modules or densely connected groups?
-   Are TLR, NF-κB or JAK-STAT-related genes prominent in the network?
-   Are there isolated nodes or small components?
-   Oriented edges ?


### 1.2 Non-inflamed samples

![PC network inferred from non-inflamed
samples](results/figures/PC_network_noninflamed.png)

| gene     | degree  | 
| -------- | ------- |
| BTK      | 6       |
| TLR8     | 6       |
| LY86     | 6       |
| CTSK     | 6       |
| LY96     | 6       |


**Interpretation**
-   Are there visibly distinct modules or densely connected groups?
-   Are TLR, NF-κB or JAK-STAT-related genes prominent in the network?
-   Are there isolated nodes or small components?
-   Oriented edges ?

### 1.3 Initial visual comparison

![Comparison of the inflamed and non-inflamed
networks](results/analysis/figures/network_comparison.png)

**Points to discuss**

-   Are the two networks globally similar, or do they show visibly
    different structures?
-   Which genes or groups of genes appear to be connected differently?
-   Are any differences concentrated around inflammation-related genes?
-   Does the figure show shared and condition-specific edges clearly?





------------------------------------------------------------------------

## 2. Global network properties

The global statistics are stored in:

-   [`network_global_statistics.csv`](analysis/tables/network_global_statistics.csv)
-   [`network_similarity.csv`](analysis/tables/network_similarity.csv)
-   [`edge_comparison_summary.csv`](analysis/tables/edge_comparison_summary.csv)

The cleaned edge lists are also available:

-   [`edges_inflamed_clean.csv`](analysis/tables/edges_inflamed_clean.csv)
-   [`edges_noninflamed_clean.csv`](analysis/tables/edges_noninflamed_clean.csv)
-   [`shared_edges.csv`](analysis/tables/shared_edges.csv)
-   [`inflamed_specific_edges.csv`](analysis/tables/inflamed_specific_edges.csv)
-   [`noninflamed_specific_edges.csv`](analysis/tables/noninflamed_specific_edges.csv)

### 2.1 Network-level statistics

Open `network_global_statistics.csv` and compare the reported metrics
for the two conditions.

  -----------------------------------------------------------------------
  Metric                  What it describes       Interpretation to check
  ----------------------- ----------------------- -----------------------
  Number of nodes         Genes represented in    Are the node sets
                          the graph               identical, as expected
                                                  from the shared
                                                  selected gene set?
                                                  Check whether isolated
                                                  nodes are retained.

  Number of edges         Inferred relationships  Which condition has
                                                  more inferred edges?

  Density                 Proportion of possible  Is one graph more
                          edges that are present  densely connected?

  Number of connected     Separate connected      Is the network more
  components              subgraphs               fragmented in one
                                                  condition?

  Largest component       Size of the largest     Is most of the network
                          connected subgraph      concentrated in one
                                                  component?

  Degree-related          Number of edges         Is connectivity
  statistics              incident to nodes       concentrated in a few
                                                  nodes or more evenly
                                                  distributed?
  -----------------------------------------------------------------------


**Results**

-   Inflamed network: 
-   Non-inflamed network: 
-   Main difference: 


### 2.2 Network similarity and edge overlap

Consult `network_similarity.csv` and `edge_comparison_summary.csv`.

-   Number of shared edges:
-   Number of inflamed-specific edges: 
-   Number of non-inflamed-specific edges:
-   Similarity/overlap metric and its definition: 



**Interpretation**



------------------------------------------------------------------------

## 3. Node-level connectivity

The node degree tables are:

-   [`node_degree_inflamed.csv`](analysis/tables/node_degree_inflamed.csv)
-   [`node_degree_noninflamed.csv`](analysis/tables/node_degree_noninflamed.csv)
-   [`node_degree_comparison.csv`](analysis/tables/node_degree_comparison.csv)

### 3.1 Highly connected genes


-   Top genes in the inflamed network: 
-   Top genes in the non-inflamed network: 
-   Genes with the largest degree differences: 

| gene     | degree inflamed  | degree non inflamed  | degree difference | 
| -------- | -------          | -------              | -------           |
| IFNG     | 5                | 2                    | 3                 |
| BTK      | 8                | 6                    | 2                 |
| LCP1     | 7                | 5                    | 2                 |
| IL10     | 7                | 5                    | 2                 |
| DUSP4    | 6                | 4                    | 2                 |
| NOD2     | 6                | 4                    | 2                 |
| DNM1     | 4                | 2                    | 2                 |
| LY86     | 4                | 6                    | -2                |
| MSN      | 4                | 2                    | 2                 |
| BIRC3    | 3                | 5                    | -2                |


**Interpretation**



### 3.2 Genes with condition-dependent connectivity

Consult `node_degree_comparison.csv`: identify genes whose degree differs
most between the two networks.



------------------------------------------------------------------------

## 4. Condition-specific interactions

The following files contain the edge-level comparison:

-   [`inflamed_specific_edges.csv`](analysis/tables/inflamed_specific_edges.csv)
-   [`noninflamed_specific_edges.csv`](analysis/tables/noninflamed_specific_edges.csv)
-   [`shared_edges.csv`](analysis/tables/shared_edges.csv)

### 4.1 Edges specific to inflamed samples



### 4.2 Edges specific to non-inflamed samples



### 4.3 Shared edges




------------------------------------------------------------------------

## 5. Reactome pathway enrichment

The enrichment results are stored in:

-   [`enrichment_inflamed.csv`](analysis/tables/enrichment_inflamed.csv)
-   [`enrichment_noninflamed.csv`](analysis/tables/enrichment_noninflamed.csv)

The corresponding R objects are:

-   [`enrichment_inflamed.rds`](analysis/networks/enrichment_inflamed.rds)
-   [`enrichment_noninflamed.rds`](analysis/networks/enrichment_noninflamed.rds)



### 5.1 Inflamed network

Open `enrichment_inflamed.csv` and inspect the enriched Reactome terms,
gene counts and adjusted p-values (or the corresponding
multiple-testing-corrected column used by the script).

-   Number of reported terms:
-   Most significant terms: 
-   Genes contributing to notable terms: 


### 5.2 Non-inflamed network

Repeat the same inspection for `enrichment_noninflamed.csv`.

-   Number of reported terms: 
-   Most significant terms: 
-   Genes contributing to notable terms: 



### 5.3 Comparison of enrichment



------------------------------------------------------------------------

## 6. Overall interpretation


### Main observations

1.  **Network structure:** 
2.  **Connectivity:** 
3.  **Condition-specific edges:** 
4.  **Pathway enrichment:** 

### Biological interpretation



### Limitations

-   **Observational data:** inferred conditional dependencies do not
    establish causal regulation.
-   **CPDAG interpretation:** the PC algorithm returns a partially
    directed graph; edge orientations may be unresolved and should not
    be overinterpreted.
-   **Patient-level dependence:** multiple biopsies may come from the
    same patient. If this dependence is not modelled, the assumption
    that all samples are independent may be violated.
-   **Separate network estimation:** condition-specific edges are
    differences between estimated graphs, not automatically
    statistically significant differential interactions.
-   **Gene selection:** both networks use a preselected set enriched for
    inflammation-related genes, so conclusions apply to this selected
    set rather than the whole transcriptome.
-   **Enrichment bias:** enrichment results are exploratory because
    candidate genes were selected partly using the same pathway
    families.
-   **Estimation stability:** the inferred edges may be sensitive to
    sample composition, sample size, preprocessing choices and PC
    algorithm parameters. Stability analysis or resampling would help
    assess robustness.

------------------------------------------------------------------------

## 7. Reproducibility and result files

All paths below are relative to the project root.

### Inferred networks

-   `results/networks/PC_network_inflamed.rds`
-   `results/networks/PC_network_noninflamed.rds`
-   `results/networks/PC_network_inflamed_edges.csv`
-   `results/networks/PC_network_noninflamed_edges.csv`

### Network figures

-   `results/figures/PC_network_inflamed.png`
-   `results/figures/PC_network_noninflamed.png`

### Analysis figures

-   `results/analysis/figures/network_comparison.png`
-   `results/analysis/figures/network_inflamed.png`
-   `results/analysis/figures/network_noninflamed.png`

### Analysis tables and objects

-   `results/analysis/tables/` --- network statistics, node degrees,
    edge comparisons and Reactome enrichment tables.
-   `results/analysis/networks/` --- RDS objects containing the analysed
    networks and enrichment results.


