# =============================================================================
# IBD GENE NETWORKS: GENE SELECTION, PC INFERENCE AND VISUALISATION
# RESYS - M2 Bioinformatique et Modélisation, 2026-2027
# =============================================================================
#
# PURPOSE
#   1. Select genes from relevant Reactome pathways.
#   2. Retain all transcription factors among those pathway genes, plus the
#      100 most variable pathway genes (variability corrected for tissue region).
#   3. Infer separate PC networks for inflamed and non-inflamed biopsies.
#   4. Export the edge lists, fitted PC objects, expression matrices and plots.
#
# INPUT
#   data/preprocessed_data.rds
#
# OUTPUT
#   data/selected_genes_annotation_PC.csv
#   data/selected_expression_matrices.rds
#   data/PC_network_inflamed_edges.csv
#   data/PC_network_noninflamed_edges.csv
#   data/PC_network_inflamed.rds
#   data/PC_network_noninflamed.rds
#   data/PC_network_inflamed.png
#   data/PC_network_noninflamed.png
#
# IMPORTANT
#   - This script assumes that preprocessing has already been completed and
#     that data/preprocessed_data.rds exists.
#   - Expression matrices have samples in rows and Ensembl gene IDs in columns.
#   - Region-corrected expression is used only to rank genes by variability.
#     Original log-CPM values are used for PC network inference.
#   - PC estimates a CPDAG (a partially directed graph). Edge orientations
#     should not be interpreted as proof of causal regulation.
#   - Samples from the same patient may not be independent. This is a
#     methodological limitation to consider when interpreting the networks.
# =============================================================================


# =============================================================================
# 0. INSTALL AND LOAD REQUIRED PACKAGES
# =============================================================================
#
# Run this section once on a new R installation. Packages already installed
# will not be reinstalled. Bioconductor packages are installed through
# BiocManager; the remaining packages are installed from CRAN.
#
# Note: installation requires an internet connection and sufficient disk space.

cran_packages = c(
  "data.table",
  "dplyr",
  "matrixStats",
  "msigdbr",
  "igraph",
  "ggraph",
  "ggplot2"
)

bioc_packages = c(
  "pcalg",
  "dorothea"
)

# Install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Use the standard Bioconductor repositories
options(repos = BiocManager::repositories())

# Install missing CRAN packages
missing_cran = cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# Install missing Bioconductor packages
missing_bioc = bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

# Load packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(matrixStats)
  library(msigdbr)
  library(pcalg)
  library(dorothea)
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

# Reproducibility
set.seed(123)


# =============================================================================
# 1. CHECK INPUT FILES AND LOAD PREPROCESSED DATA
# =============================================================================

input_file = "data/processed/preprocessed_data.rds"

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nRun the preprocessing script first."
  )
}

preprocessed = readRDS(input_file)

# Extract the three objects saved by the preprocessing script
expr_clean = as.data.table(preprocessed$expr_clean)
annotation_kept = as.data.table(preprocessed$annotation_kept)
final_samples = as.data.table(preprocessed$final_samples)

cat("\nINPUT DATA\n")
cat("Expression table:", nrow(expr_clean), "samples x",
    ncol(expr_clean) - 1, "genes\n")
cat("Gene annotation:", nrow(annotation_kept), "rows\n")
cat("Sample metadata:", nrow(final_samples), "rows\n")

# Check that the expected columns are available
if (!"samples_name" %in% names(expr_clean)) {
  stop("Column 'samples_name' is missing from expr_clean.")
}

if (!all(c("samples_name", "typere", "region_groups") %in%
         names(final_samples))) {
  stop(
    "final_samples must contain: samples_name, typere and region_groups."
  )
}

if (!all(c("EnsemblGeneID", "Symbol") %in% names(annotation_kept))) {
  stop(
    "annotation_kept must contain: EnsemblGeneID and Symbol."
  )
}

# Convert expression values to a numeric matrix.
# Rows = samples; columns = Ensembl gene IDs.
expr = as.data.frame(
  expr_clean[, setdiff(names(expr_clean), "samples_name"), with = FALSE]
)

rownames(expr) = expr_clean$samples_name
expr = as.matrix(expr)
storage.mode(expr) = "numeric"

# Align metadata rows with the expression matrix.
sample_index = match(rownames(expr), final_samples$samples_name)

if (anyNA(sample_index)) {
  stop("Some expression samples are missing from final_samples.")
}

final_samples = final_samples[sample_index]

stopifnot(
  identical(rownames(expr), final_samples$samples_name),
  !anyDuplicated(rownames(expr)),
  !anyDuplicated(colnames(expr))
)

if (anyNA(expr) || any(!is.finite(expr))) {
  stop("The expression matrix contains missing or non-finite values.")
}

cat("\nAligned expression matrix:", nrow(expr), "samples x",
    ncol(expr), "genes\n")

cat("\nSamples by inflammatory status:\n")
print(table(final_samples$typere, useNA = "ifany"))

cat("\nSamples by anatomical region:\n")
print(table(final_samples$region_groups, useNA = "ifany"))


# =============================================================================
# 2. BIOLOGICAL GENE SELECTION
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Retrieve Reactome pathways from MSigDB
# -----------------------------------------------------------------------------

# msigdbr changed its argument names in recent versions. This block supports
# both the current API (collection/subcollection) and the older API
# (category/subcategory).

msig_formals = names(formals(msigdbr::msigdbr))

if ("collection" %in% msig_formals) {
  reactome = msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"
  )
} else {
  reactome = msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "C2",
    subcategory = "CP:REACTOME"
  )
}

# Select pathways whose names match the biological themes of interest.
# Inspect the printed list: pathway names in MSigDB may change between releases.
pathway_pattern = paste(
  "TOLL_LIKE_RECEPTOR",
  "NFKB",
  "JAK_STAT",
  sep = "|"
)

selected_pathways = reactome %>%
  filter(grepl(pathway_pattern, gs_name, ignore.case = TRUE))

selected_pathway_names = sort(unique(selected_pathways$gs_name))

cat("\nSELECTED REACTOME PATHWAYS\n")
print(selected_pathway_names)

if (nrow(selected_pathways) == 0) {
  stop(
    "No Reactome pathways matched the requested pattern. ",
    "Inspect the available pathway names in the current msigdbr release."
  )
}

# Unique gene symbols belonging to the selected pathways
pathway_symbols = unique(selected_pathways$gene_symbol)
pathway_symbols = pathway_symbols[!is.na(pathway_symbols)]

cat("\nUnique symbols in selected pathways:",
    length(pathway_symbols), "\n")


# -----------------------------------------------------------------------------
# 2.2 Retrieve human transcription factors from DoRothEA
# -----------------------------------------------------------------------------

# DoRothEA provides a curated catalogue of transcription factors.
# Here, we retain TFs only if they also belong to the selected pathways.

data("dorothea_hs", package = "dorothea", envir = environment())

if (!exists("dorothea_hs")) {
  stop("The dataset 'dorothea_hs' could not be loaded from dorothea.")
}

tf_symbols = unique(dorothea_hs$tf)
tf_symbols = tf_symbols[!is.na(tf_symbols)]

cat("Human TF symbols in DoRothEA:", length(tf_symbols), "\n")


# -----------------------------------------------------------------------------
# 2.3 Map gene symbols to Ensembl IDs
# -----------------------------------------------------------------------------

# Expression columns are Ensembl IDs, whereas MSigDB and DoRothEA use symbols.
# The annotation table provides the mapping between these identifiers.

gene_annotation = unique(
  annotation_kept[
    !is.na(EnsemblGeneID) &
      EnsemblGeneID != "" &
      !is.na(Symbol) &
      Symbol != "",
    .(EnsemblGeneID, Symbol)
  ]
)

# Pathway genes that are actually present in the expression matrix
pathway_genes = gene_annotation[
  Symbol %in% pathway_symbols &
    EnsemblGeneID %in% colnames(expr),
  unique(EnsemblGeneID)
]

# Transcription factors among those pathway genes
pathway_tfs = gene_annotation[
  Symbol %in% tf_symbols &
    EnsemblGeneID %in% pathway_genes,
  unique(EnsemblGeneID)
]

cat("\nPathway genes present in expression matrix:",
    length(pathway_genes), "\n")
cat("TFs among selected pathway genes:",
    length(pathway_tfs), "\n")

if (length(pathway_genes) == 0) {
  stop("No selected pathway genes were found in the expression matrix.")
}


# =============================================================================
# 3. REGION-CORRECTED VARIABILITY (MEDIAN ABSOLUTE DEVIATION)
# =============================================================================

# Anatomical region is a potential source of expression variation.
# We remove its linear effect before ranking pathway genes by variability.
#
# IMPORTANT:
#   - This correction is used only for gene selection.
#   - The original log-CPM values are retained for PC inference.
#   - This is a linear adjustment for region; it does not remove every
#     possible confounding factor.

region = factor(final_samples$region_groups)

if (anyNA(region)) {
  stop("Missing values were found in region_groups.")
}

# Design matrix: intercept + anatomical region
design_region = model.matrix(~ region)

# Fit a linear model separately for each gene and retain its residuals.
expr_region_corrected = vapply(
  seq_len(ncol(expr)),
  function(j) {
    lm.fit(
      x = design_region,
      y = expr[, j]
    )$residuals
  },
  numeric(nrow(expr))
)

rownames(expr_region_corrected) = rownames(expr)
colnames(expr_region_corrected) = colnames(expr)

# Calculate the median absolute deviation (MAD) for each gene
gene_mad = matrixStats::colMads(
  expr_region_corrected,
  na.rm = TRUE
)

mad_table = data.table(
  EnsemblGeneID = names(gene_mad),
  MAD = as.numeric(gene_mad)
)

# Rank only pathway genes by their region-corrected MAD
mad_pathway = mad_table[
  EnsemblGeneID %in% pathway_genes
]

setorder(mad_pathway, -MAD)

cat("\nTOP 10 PATHWAY GENES BY REGION-CORRECTED MAD\n")
print(head(mad_pathway, 10))


# =============================================================================
# 4. FINAL GENE SELECTION
# =============================================================================

# Number of most variable pathway genes to retain.
# This is an analysis parameter and can be adjusted after sensitivity checks.
top_mad_n = 100

top_variable_genes = head(
  mad_pathway$EnsemblGeneID,
  top_mad_n
)

# Final selection:
#   1. All transcription factors belonging to the selected pathways.
#   2. The top 100 most variable genes among the selected pathway genes.

selected_genes = unique(c(
  pathway_tfs,
  top_variable_genes
))

# Keep genes present in the expression matrix
selected_genes = intersect(selected_genes, colnames(expr))

if (length(selected_genes) < 2) {
  stop("Fewer than two genes were selected; PC inference cannot proceed.")
}

# Keep annotation rows in the same order as selected_genes.
selected_annotation = gene_annotation[
  match(selected_genes, EnsemblGeneID)
]

# Add MAD values without changing the selected gene order
selected_annotation[
  mad_table,
  MAD := i.MAD,
  on = "EnsemblGeneID"
]

cat("\nFINAL GENE SELECTION\n")
cat("Pathway genes available:", length(pathway_genes), "\n")
cat("Pathway transcription factors:", length(pathway_tfs), "\n")
cat("Top variable pathway genes:", length(top_variable_genes), "\n")
cat("Unique genes selected:", length(selected_genes), "\n")

fwrite(
  selected_annotation,
  "data/processed/selected_genes_annotation_PC.csv"
)


# =============================================================================
# 5. BUILD THE TWO EXPRESSION SUBMATRICES
# =============================================================================

# The same selected genes are used in both conditions so that the two
# inferred networks are directly comparable.

# Select columns in the expression matrix
expr_selected = expr[, selected_genes, drop = FALSE]

# Identify samples by inflammatory status
inflamed_samples = final_samples$samples_name[
  final_samples$typere == "I"
]

noninflamed_samples = final_samples$samples_name[
  final_samples$typere == "NonI"
]

if (length(inflamed_samples) < 3 || length(noninflamed_samples) < 3) {
  stop("Too few samples in at least one inflammatory-status group.")
}

# Subset expression matrices
expr_inflamed = expr_selected[
  inflamed_samples,
  ,
  drop = FALSE
]

expr_noninflamed = expr_selected[
  noninflamed_samples,
  ,
  drop = FALSE
]

cat("\nInflamed expression matrix:",
    nrow(expr_inflamed), "samples x", ncol(expr_inflamed), "genes\n")
cat("Non-inflamed expression matrix:",
    nrow(expr_noninflamed), "samples x", ncol(expr_noninflamed), "genes\n")

# Retain genes with non-zero variance in BOTH groups.
# This ensures that both networks are inferred on exactly the same variables.
variance_inflamed = apply(expr_inflamed, 2, var)
variance_noninflamed = apply(expr_noninflamed, 2, var)

valid_genes = names(variance_inflamed)[
  is.finite(variance_inflamed) &
    is.finite(variance_noninflamed) &
    variance_inflamed > 0 &
    variance_noninflamed > 0
]

expr_inflamed = expr_inflamed[, valid_genes, drop = FALSE]
expr_noninflamed = expr_noninflamed[, valid_genes, drop = FALSE]

if (length(valid_genes) < 2) {
  stop("Fewer than two variable genes remain in both conditions.")
}

stopifnot(identical(colnames(expr_inflamed), colnames(expr_noninflamed)))

cat("\nGenes retained for both networks:", length(valid_genes), "\n")

# Keep the annotation restricted to the final set of genes
selected_annotation = selected_annotation[
  match(valid_genes, EnsemblGeneID)
]

# Save matrices and their annotation for subsequent scripts
saveRDS(
  list(
    expr_inflamed = expr_inflamed,
    expr_noninflamed = expr_noninflamed,
    selected_genes = valid_genes,
    selected_annotation = selected_annotation
  ),
  "data/processed/selected_expression_matrices.rds",
  compress = "gzip"
)


# =============================================================================
# 6. PC NETWORK INFERENCE
# =============================================================================

# The PC algorithm estimates a CPDAG from conditional-independence tests.
# Here we use:
#   - Gaussian conditional-independence tests (gaussCItest)
#   - alpha = 0.01
#   - m.max = 2 (conditioning sets of at most two variables)
#
# These are analysis choices, not universal defaults. Their influence should
# ideally be assessed in a sensitivity analysis.

infer_PC = function(expression, alpha = 0.01, m.max = 2) {
  
  expression = as.matrix(expression)
  
  n_samples = nrow(expression)
  n_genes = ncol(expression)
  
  cat("\nInferring PC network\n")
  cat("Samples:", n_samples, "\n")
  cat("Genes:", n_genes, "\n")
  cat("alpha:", alpha, "\n")
  cat("Maximum conditioning-set size:", m.max, "\n")
  
  if (any(!is.finite(expression))) {
    stop("Expression matrix contains non-finite values.")
  }
  
  # Pearson correlation matrix across genes
  correlation = cor(
    expression,
    method = "pearson"
  )
  
  if (any(!is.finite(correlation))) {
    stop("Non-finite values were found in the correlation matrix.")
  }
  
  # Run the PC algorithm
  pc_fit = pcalg::pc(
    suffStat = list(
      C = correlation,
      n = n_samples
    ),
    indepTest = pcalg::gaussCItest,
    labels = colnames(expression),
    alpha = alpha,
    skel.method = "stable.fast",
    m.max = m.max
  )
  
  return(pc_fit)
}


# -----------------------------------------------------------------------------
# 6.1 Infer the network in inflamed tissue
# -----------------------------------------------------------------------------

pc_inflamed = infer_PC(
  expression = expr_inflamed,
  alpha = 0.01,
  m.max = 2
)


# -----------------------------------------------------------------------------
# 6.2 Infer the network in non-inflamed tissue
# -----------------------------------------------------------------------------

pc_noninflamed = infer_PC(
  expression = expr_noninflamed,
  alpha = 0.01,
  m.max = 2
)


# =============================================================================
# 7. EXTRACT AND EXPORT NETWORK EDGES
# =============================================================================

# The PC result is a partially directed graph (CPDAG).
# We retain the edge type in the exported table:
#   - "undirected": an undirected edge in the CPDAG
#   - "directed": an oriented edge in the CPDAG
#
# The 'from' and 'to' columns are identifiers used to encode the edge.
# For undirected edges, their order has no biological meaning.

extract_edges = function(pc_fit) {
  
  adjacency = as(
    pc_fit@graph,
    "matrix"
  )
  
  genes = colnames(adjacency)
  
  # Each pair is extracted once, using the upper triangle.
  edge_indices = which(
    upper.tri(adjacency) &
      (
        adjacency != 0 |
          t(adjacency) != 0
      ),
    arr.ind = TRUE
  )
  
  if (nrow(edge_indices) == 0) {
    return(data.table(
      from = character(),
      to = character(),
      type = character()
    ))
  }
  
  edges = data.table(
    from = genes[edge_indices[, 1]],
    to = genes[edge_indices[, 2]],
    type = vapply(
      seq_len(nrow(edge_indices)),
      function(i) {
        
        row_i = edge_indices[i, 1]
        col_i = edge_indices[i, 2]
        
        a = adjacency[row_i, col_i]
        b = adjacency[col_i, row_i]
        
        if (a != 0 && b != 0) {
          "undirected"
        } else {
          "directed"
        }
      },
      character(1)
    )
  )
  
  return(edges)
}

edges_inflamed = extract_edges(pc_inflamed)
edges_noninflamed = extract_edges(pc_noninflamed)

# Export edge lists
fwrite(
  edges_inflamed,
  "results/networks/PC_network_inflamed_edges.csv"
)

fwrite(
  edges_noninflamed,
  "results/networks/PC_network_noninflamed_edges.csv"
)

# Save fitted PC objects for later inspection
saveRDS(
  pc_inflamed,
  "results/networks/PC_network_inflamed.rds"
)

saveRDS(
  pc_noninflamed,
  "results/networks/PC_network_noninflamed.rds"
)


# =============================================================================
# 8. BUILD IG RAPH OBJECTS FOR VISUALISATION
# =============================================================================

# Keep Ensembl IDs as graph vertex names. This avoids accidentally merging
# distinct Ensembl IDs that happen to share the same gene symbol.
#
# Include every selected gene as a vertex, including genes with no inferred
# edges. Such genes are isolates and would otherwise disappear from the graph.

gene_labels = selected_annotation$Symbol
names(gene_labels) = selected_annotation$EnsemblGeneID

vertices = data.frame(
  name = valid_genes,
  symbol = unname(gene_labels[valid_genes]),
  stringsAsFactors = FALSE
)

# Replace missing symbols with Ensembl IDs for readable plot labels
missing_symbols = is.na(vertices$symbol) | vertices$symbol == ""
vertices$symbol[missing_symbols] = vertices$name[missing_symbols]

network_inflamed = igraph::graph_from_data_frame(
  d = edges_inflamed[, .(from, to, type)],
  directed = FALSE,
  vertices = vertices
)

network_noninflamed = igraph::graph_from_data_frame(
  d = edges_noninflamed[, .(from, to, type)],
  directed = FALSE,
  vertices = vertices
)

cat("\nNETWORK OBJECTS\n")
cat("Inflamed: ", igraph::vcount(network_inflamed), "vertices, ",
    igraph::ecount(network_inflamed), "edges\n", sep = "")
cat("Non-inflamed: ", igraph::vcount(network_noninflamed), "vertices, ",
    igraph::ecount(network_noninflamed), "edges\n", sep = "")


# =============================================================================
# 9. VISUALISE THE INFLAMED NETWORK
# =============================================================================

set.seed(123)

p_inflamed = ggraph(
  network_inflamed,
  layout = "fr"
) +
  geom_edge_link(
    colour = "grey75",
    alpha = 0.6,
    linewidth = 0.5
  ) +
  geom_node_point(
    colour = "#B65C54",
    size = 3
  ) +
  geom_node_text(
    aes(label = symbol),
    repel = TRUE,
    size = 3
  ) +
  labs(
    title = "Gene network - Inflamed tissue"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    )
  )

print(p_inflamed)

ggsave(
  filename = "results/figures/PC_network_inflamed.png",
  plot = p_inflamed,
  width = 12,
  height = 10,
  dpi = 300
)


# =============================================================================
# 10. VISUALISE THE NON-INFLAMED NETWORK
# =============================================================================

set.seed(123)

p_noninflamed = ggraph(
  network_noninflamed,
  layout = "fr"
) +
  geom_edge_link(
    colour = "grey75",
    alpha = 0.6,
    linewidth = 0.5
  ) +
  geom_node_point(
    colour = "#3B7896",
    size = 3
  ) +
  geom_node_text(
    aes(label = symbol),
    repel = TRUE,
    size = 3
  ) +
  labs(
    title = "Gene network - Non-inflamed tissue"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    )
  )

print(p_noninflamed)

ggsave(
  filename = "results/figures/PC_network_noninflamed.png",
  plot = p_noninflamed,
  width = 12,
  height = 10,
  dpi = 300
)


# =============================================================================
# 11. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("PC NETWORK INFERENCE COMPLETED\n")
cat("============================================================\n")

cat("\nGene selection\n")
cat("  Pathway genes available:", length(pathway_genes), "\n")
cat("  Pathway TFs:", length(pathway_tfs), "\n")
cat("  Final genes used in both networks:", length(valid_genes), "\n")

cat("\nInflamed network\n")
cat("  Samples:", nrow(expr_inflamed), "\n")
cat("  Vertices:", igraph::vcount(network_inflamed), "\n")
cat("  Edges:", igraph::ecount(network_inflamed), "\n")

cat("\nNon-inflamed network\n")
cat("  Samples:", nrow(expr_noninflamed), "\n")
cat("  Vertices:", igraph::vcount(network_noninflamed), "\n")
cat("  Edges:", igraph::ecount(network_noninflamed), "\n")

cat("\nFiles saved in data/\n")
cat("  - selected_genes_annotation_PC.csv\n")
cat("  - selected_expression_matrices.rds\n")
cat("  - PC_network_inflamed_edges.csv\n")
cat("  - PC_network_noninflamed_edges.csv\n")
cat("  - PC_network_inflamed.rds\n")
cat("  - PC_network_noninflamed.rds\n")
cat("  - PC_network_inflamed.png\n")
cat("  - PC_network_noninflamed.png\n")

cat("\nReminder: PC edge orientations do not establish causality.\n")
cat("============================================================\n")
