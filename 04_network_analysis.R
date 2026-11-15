# ============================================================
# COMPARISON AND ANALYSIS OF INFLAMED / NON-INFLAMED NETWORKS
# ============================================================
#
# Objectif :
#   - comparer les réseaux PC inférés dans les tissus inflammés
#     et non inflammés ;
#   - décrire leur topologie et les différences d'interactions ;
#   - visualiser les réseaux et exporter les résultats ;
#   - réaliser un enrichissement Reactome exploratoire.
#
# Prérequis :
#   - avoir exécuté le script d'inférence des réseaux PC ;
#   - travailler depuis la racine du projet ;
#   - disposer des fichiers dans le dossier data/.
#
# Tous les résultats de ce script sont enregistrés dans :
#   results/network_analysis/
#
# IMPORTANT : une arête détectée dans un seul réseau n'est pas,
# à elle seule, la preuve d'une différence statistiquement
# significative entre les deux conditions.
# ============================================================


# ============================================================
# 0. INSTALLATION AND LOADING OF DEPENDENCIES
# ============================================================

# Packages CRAN
cran_packages = c(
  "data.table",
  "dplyr",
  "igraph",
  "ggraph",
  "ggplot2",
  "ggrepel",
  "msigdbr"
)

# Packages Bioconductor
bioc_packages = c(
  "clusterProfiler",
  "GO.db"
)

# Install only missing CRAN packages
missing_cran = cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# Configure Bioconductor repositories
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

options(repos = BiocManager::repositories())

# Install only missing Bioconductor packages
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
library(data.table)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(clusterProfiler)


# ============================================================
# 1. PROJECT PATHS AND OUTPUT FOLDERS
# ============================================================

# Input files
input_edges_inflamed = "results/networks/PC_network_inflamed_edges.csv"
input_edges_noninflamed = "results/networks/PC_network_noninflamed_edges.csv"
input_selected_data = "data/processed/selected_expression_matrices.rds"

# Main output directory
output_dir = "results/analysis"

# Subdirectories
output_tables = file.path(output_dir, "tables")
output_figures = file.path(output_dir, "figures")
output_networks = file.path(output_dir, "networks")

dir.create(output_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(output_networks, recursive = TRUE, showWarnings = FALSE)

# Check that all input files exist
input_files = c(
  input_edges_inflamed,
  input_edges_noninflamed,
  input_selected_data
)

missing_files = input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "The following input files are missing:\n",
    paste(missing_files, collapse = "\n"),
    "\nRun the network inference script first, and check the project directory."
  )
}


# ============================================================
# 2. LOAD DATA AND GENE ANNOTATIONS
# ============================================================

edges_inflamed = fread(input_edges_inflamed)
edges_noninflamed = fread(input_edges_noninflamed)

selected_data = readRDS(input_selected_data)

if (is.null(selected_data$selected_annotation)) {
  stop("selected_expression_matrices.rds must contain 'selected_annotation'.")
}

selected_annotation = as.data.table(selected_data$selected_annotation)

required_annotation_columns = c("EnsemblGeneID", "Symbol")

if (!all(required_annotation_columns %in% names(selected_annotation))) {
  stop(
    "selected_annotation must contain the columns: ",
    paste(required_annotation_columns, collapse = ", ")
  )
}

# Ensembl ID -> gene symbol mapping
gene_labels = selected_annotation$Symbol
names(gene_labels) = selected_annotation$EnsemblGeneID

# Remove duplicated identifiers from the mapping, if any
gene_labels = gene_labels[!duplicated(names(gene_labels))]


# ============================================================
# 3. PREPARE EDGE LISTS
# ============================================================

# The PC networks are compared as undirected networks.
# Thus, A-B and B-A are treated as the same interaction.

prepare_edges = function(edges, condition_name) {
  
  edges = as.data.table(edges)
  
  if (!all(c("from", "to") %in% names(edges))) {
    stop(
      "The edge table for ", condition_name,
      " must contain columns named 'from' and 'to'."
    )
  }
  
  edges[, from := as.character(from)]
  edges[, to := as.character(to)]
  
  # Remove missing endpoints and self-loops
  edges = edges[
    !is.na(from) & !is.na(to) & from != to
  ]
  
  # Canonical identifier for undirected edges
  edges[, edge_id := paste(
    pmin(from, to),
    pmax(from, to),
    sep = "__"
  )]
  
  # Keep one row per undirected edge
  edges = unique(edges, by = "edge_id")
  
  # Add gene symbols for interpretation
  edges[, from_symbol := unname(gene_labels[from])]
  edges[, to_symbol := unname(gene_labels[to])]
  
  # Keep a condition label for downstream exports
  edges[, condition := condition_name]
  
  return(edges)
}

edges_inflamed = prepare_edges(edges_inflamed, "Inflamed")
edges_noninflamed = prepare_edges(edges_noninflamed, "Non-inflamed")

cat("\nNumber of edges after cleaning:\n")
cat("Inflamed:", nrow(edges_inflamed), "\n")
cat("Non-inflamed:", nrow(edges_noninflamed), "\n")

# Export cleaned edge lists
fwrite(
  edges_inflamed,
  file.path(output_tables, "edges_inflamed_clean.csv")
)

fwrite(
  edges_noninflamed,
  file.path(output_tables, "edges_noninflamed_clean.csv")
)


# ============================================================
# 4. BUILD NETWORK OBJECTS
# ============================================================

# Only genes participating in at least one edge are included here.
# Isolated genes are therefore not represented in these graphs.

network_inflamed = graph_from_data_frame(
  edges_inflamed[, .(from, to)],
  directed = FALSE
)

network_noninflamed = graph_from_data_frame(
  edges_noninflamed[, .(from, to)],
  directed = FALSE
)

# Store gene symbols as vertex attributes.
V(network_inflamed)$symbol = unname(gene_labels[V(network_inflamed)$name])
V(network_noninflamed)$symbol = unname(gene_labels[V(network_noninflamed)$name])

# Save graph objects so they can be reused without rebuilding them.
saveRDS(
  network_inflamed,
  file.path(output_networks, "network_inflamed.rds")
)

saveRDS(
  network_noninflamed,
  file.path(output_networks, "network_noninflamed.rds")
)


# ============================================================
# 5. GLOBAL NETWORK STATISTICS
# ============================================================

network_stats = function(network, condition_name) {
  
  n_nodes = vcount(network)
  n_edges = ecount(network)
  
  # Some metrics are undefined for an empty graph.
  if (n_nodes == 0) {
    return(data.frame(
      condition = condition_name,
      n_nodes = 0,
      n_edges = 0,
      density = NA_real_,
      n_components = 0,
      mean_degree = NA_real_,
      max_degree = NA_real_,
      mean_clustering = NA_real_
    ))
  }
  
  degrees = degree(network)
  
  data.frame(
    condition = condition_name,
    n_nodes = n_nodes,
    n_edges = n_edges,
    density = edge_density(network),
    n_components = igraph::components(network)$no,
    mean_degree = mean(degrees),
    max_degree = max(degrees),
    mean_clustering = transitivity(
      network,
      type = "average",
      isolates = "zero"
    )
  )
}

network_global_statistics = bind_rows(
  network_stats(network_inflamed, "Inflamed"),
  network_stats(network_noninflamed, "Non-inflamed")
)

print(network_global_statistics)

fwrite(
  network_global_statistics,
  file.path(output_tables, "network_global_statistics.csv")
)


# ============================================================
# 6. NODE DEGREE AND HUB GENES
# ============================================================

# Degree = number of edges connected to a node.
# It is a simple measure of connectivity, not a measure of
# biological importance or regulatory activity.

get_node_degree = function(network, condition_name) {
  
  if (vcount(network) == 0) {
    return(data.frame(
      gene_id = character(),
      symbol = character(),
      degree = integer(),
      condition = character()
    ))
  }
  
  data.frame(
    gene_id = V(network)$name,
    symbol = V(network)$symbol,
    degree = as.integer(degree(network)),
    condition = condition_name
  ) %>%
    arrange(desc(degree), gene_id)
}

degree_inflamed = get_node_degree(network_inflamed, "Inflamed")
degree_noninflamed = get_node_degree(network_noninflamed, "Non-inflamed")

# Display the 15 most connected genes in each network
cat("\nTop 15 genes by degree - inflamed:\n")
print(head(degree_inflamed, 15))

cat("\nTop 15 genes by degree - non-inflamed:\n")
print(head(degree_noninflamed, 15))

fwrite(
  degree_inflamed,
  file.path(output_tables, "node_degree_inflamed.csv")
)

fwrite(
  degree_noninflamed,
  file.path(output_tables, "node_degree_noninflamed.csv")
)


# ============================================================
# 7. COMPARE NODE DEGREES BETWEEN CONDITIONS
# ============================================================

# Only genes present in both networks can be compared directly.
degree_comparison = inner_join(
  degree_inflamed %>%
    select(gene_id, symbol, degree_inflamed = degree),
  degree_noninflamed %>%
    select(gene_id, degree_noninflamed = degree),
  by = "gene_id"
) %>%
  mutate(
    degree_difference = degree_inflamed - degree_noninflamed
  ) %>%
  arrange(desc(abs(degree_difference)))

cat("\nGenes with the largest absolute degree differences:\n")
print(head(degree_comparison, 20))

fwrite(
  degree_comparison,
  file.path(output_tables, "node_degree_comparison.csv")
)


# ============================================================
# 8. COMPARE NETWORK EDGES
# ============================================================

shared_edges = intersect(
  edges_inflamed$edge_id,
  edges_noninflamed$edge_id
)

inflamed_specific = setdiff(
  edges_inflamed$edge_id,
  edges_noninflamed$edge_id
)

noninflamed_specific = setdiff(
  edges_noninflamed$edge_id,
  edges_inflamed$edge_id
)

edge_comparison_summary = data.frame(
  category = c(
    "Shared",
    "Inflamed-specific",
    "Non-inflamed-specific"
  ),
  n_edges = c(
    length(shared_edges),
    length(inflamed_specific),
    length(noninflamed_specific)
  )
)

print(edge_comparison_summary)

fwrite(
  edge_comparison_summary,
  file.path(output_tables, "edge_comparison_summary.csv")
)

# Jaccard similarity:
# 1 = identical edge sets; 0 = no shared edges.
all_edge_ids = union(
  edges_inflamed$edge_id,
  edges_noninflamed$edge_id
)

jaccard_edges = if (length(all_edge_ids) == 0) {
  NA_real_
} else {
  length(shared_edges) / length(all_edge_ids)
}

# Node overlap, considering only genes represented by at least
# one edge in each graph.
nodes_inflamed = V(network_inflamed)$name
nodes_noninflamed = V(network_noninflamed)$name

all_node_ids = union(nodes_inflamed, nodes_noninflamed)

jaccard_nodes = if (length(all_node_ids) == 0) {
  NA_real_
} else {
  length(intersect(nodes_inflamed, nodes_noninflamed)) /
    length(all_node_ids)
}

network_similarity = data.frame(
  metric = c("Edge Jaccard similarity", "Node Jaccard similarity"),
  value = c(jaccard_edges, jaccard_nodes)
)

print(network_similarity)

fwrite(
  network_similarity,
  file.path(output_tables, "network_similarity.csv")
)


# ============================================================
# 9. EXPORT SHARED AND CONDITION-SPECIFIC EDGES
# ============================================================

inflamed_specific_edges = edges_inflamed[
  edge_id %in% inflamed_specific
]

noninflamed_specific_edges = edges_noninflamed[
  edge_id %in% noninflamed_specific
]

shared_edges_table = edges_inflamed[
  edge_id %in% shared_edges
]

fwrite(
  inflamed_specific_edges,
  file.path(output_tables, "inflamed_specific_edges.csv")
)

fwrite(
  noninflamed_specific_edges,
  file.path(output_tables, "noninflamed_specific_edges.csv")
)

fwrite(
  shared_edges_table,
  file.path(output_tables, "shared_edges.csv")
)


# ============================================================
# 10. VISUALISE THE TWO NETWORKS
# ============================================================

# A fixed seed makes the force-directed layout reproducible.
set.seed(123)

if (vcount(network_inflamed) > 0) {
  
  p_inflamed = ggraph(network_inflamed, layout = "fr") +
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
    labs(title = "Gene network – Inflamed tissue") +
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
    file.path(output_figures, "network_inflamed.png"),
    p_inflamed,
    width = 12,
    height = 10,
    dpi = 300
  )
}

if (vcount(network_noninflamed) > 0) {
  
  p_noninflamed = ggraph(network_noninflamed, layout = "fr") +
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
    labs(title = "Gene network – Non-inflamed tissue") +
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
    file.path(output_figures, "network_noninflamed.png"),
    p_noninflamed,
    width = 12,
    height = 10,
    dpi = 300
  )
}


# ============================================================
# 11. BUILD AND VISUALISE THE COMPARISON NETWORK
# ============================================================

# Each edge is labelled according to its presence in the two
# inferred networks.

all_edges = bind_rows(
  edges_inflamed[
    edge_id %in% shared_edges,
    .(from, to, status = "Shared")
  ],
  edges_inflamed[
    edge_id %in% inflamed_specific,
    .(from, to, status = "Inflamed-specific")
  ],
  edges_noninflamed[
    edge_id %in% noninflamed_specific,
    .(from, to, status = "Non-inflamed-specific")
  ]
)

comparison_network = graph_from_data_frame(
  all_edges,
  directed = FALSE
)

if (vcount(comparison_network) > 0) {
  
  V(comparison_network)$symbol =
    unname(gene_labels[V(comparison_network)$name])
  
  saveRDS(
    comparison_network,
    file.path(output_networks, "comparison_network.rds")
  )
  
  set.seed(123)
  
  p_comparison = ggraph(comparison_network, layout = "fr") +
    geom_edge_link(
      aes(edge_colour = status),
      alpha = 0.7,
      linewidth = 0.6
    ) +
    geom_node_point(
      colour = "grey25",
      size = 3
    ) +
    geom_node_text(
      aes(label = symbol),
      repel = TRUE,
      size = 3
    ) +
    scale_edge_colour_manual(
      values = c(
        "Shared" = "grey75",
        "Inflamed-specific" = "#B65C54",
        "Non-inflamed-specific" = "#3B7896"
      )
    ) +
    labs(
      title = "Comparison of the two inferred gene networks",
      edge_colour = "Edge category"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom"
    )
  
  print(p_comparison)
  
  ggsave(
    file.path(output_figures, "network_comparison.png"),
    p_comparison,
    width = 14,
    height = 12,
    dpi = 300
  )
} else {
  warning("The comparison graph is empty: no edges were found in either network.")
}


# ============================================================
# 12. LOAD REACTOME GENE SETS
# ============================================================

# msigdbr has changed its argument names across versions.
# This block supports both the current and older interfaces.

msigdbr_arguments = names(formals(msigdbr::msigdbr))

if ("collection" %in% msigdbr_arguments) {
  
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

term2gene = reactome %>%
  select(gs_name, gene_symbol) %>%
  distinct()

if (nrow(term2gene) == 0) {
  stop("No Reactome gene sets were retrieved from msigdbr.")
}


# ============================================================
# 13. DEFINE GENE LISTS AND BACKGROUND
# ============================================================

# Enrichment is performed on genes that are endpoints of
# condition-specific edges.

genes_inflamed_specific = unique(c(
  inflamed_specific_edges$from_symbol,
  inflamed_specific_edges$to_symbol
))

genes_noninflamed_specific = unique(c(
  noninflamed_specific_edges$from_symbol,
  noninflamed_specific_edges$to_symbol
))

# Remove missing and empty symbols
genes_inflamed_specific = genes_inflamed_specific[
  !is.na(genes_inflamed_specific) & genes_inflamed_specific != ""
]

genes_noninflamed_specific = genes_noninflamed_specific[
  !is.na(genes_noninflamed_specific) & genes_noninflamed_specific != ""
]

# Background = all genes that could have entered the network
# according to the selected gene set used during inference.
gene_universe = unique(
  selected_annotation$Symbol[
    !is.na(selected_annotation$Symbol) &
      selected_annotation$Symbol != ""
  ]
)

# Restrict input gene lists to the background universe.
genes_inflamed_specific = intersect(
  genes_inflamed_specific,
  gene_universe
)

genes_noninflamed_specific = intersect(
  genes_noninflamed_specific,
  gene_universe
)

cat("\nGenes submitted to Reactome enrichment:\n")
cat("Inflamed-specific:", length(genes_inflamed_specific), "\n")
cat("Non-inflamed-specific:", length(genes_noninflamed_specific), "\n")
cat("Background:", length(gene_universe), "\n")


# ============================================================
# 14. RUN REACTOME ENRICHMENT
# ============================================================

run_reactome_enrichment = function(gene_list, universe, term2gene) {
  
  if (length(gene_list) == 0) {
    message("No genes available for enrichment in this condition.")
    return(NULL)
  }
  
  clusterProfiler::enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}

enrichment_inflamed = run_reactome_enrichment(
  genes_inflamed_specific,
  gene_universe,
  term2gene
)

enrichment_noninflamed = run_reactome_enrichment(
  genes_noninflamed_specific,
  gene_universe,
  term2gene
)

enrichment_inflamed_table = if (is.null(enrichment_inflamed)) {
  data.frame()
} else {
  as.data.frame(enrichment_inflamed)
}

enrichment_noninflamed_table = if (is.null(enrichment_noninflamed)) {
  data.frame()
} else {
  as.data.frame(enrichment_noninflamed)
}


# ============================================================
# 15. EXPORT ENRICHMENT RESULTS
# ============================================================

fwrite(
  enrichment_inflamed_table,
  file.path(output_tables, "enrichment_inflamed.csv")
)

fwrite(
  enrichment_noninflamed_table,
  file.path(output_tables, "enrichment_noninflamed.csv")
)

# Save complete enrichment objects as well.
saveRDS(
  enrichment_inflamed,
  file.path(output_networks, "enrichment_inflamed.rds")
)

saveRDS(
  enrichment_noninflamed,
  file.path(output_networks, "enrichment_noninflamed.rds")
)


# ============================================================
# 16. ENRICHMENT DOTPLOTS
# ============================================================

if (!is.null(enrichment_inflamed) &&
    nrow(enrichment_inflamed_table) > 0) {
  
  p_enrichment_inflamed = dotplot(
    enrichment_inflamed,
    showCategory = 15
  ) +
    labs(title = "Reactome enrichment – Inflamed-specific genes")
  
  print(p_enrichment_inflamed)
  
  ggsave(
    file.path(output_figures, "enrichment_inflamed.png"),
    p_enrichment_inflamed,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  message("No significant Reactome pathway found for inflamed-specific genes.")
}

if (!is.null(enrichment_noninflamed) &&
    nrow(enrichment_noninflamed_table) > 0) {
  
  p_enrichment_noninflamed = dotplot(
    enrichment_noninflamed,
    showCategory = 15
  ) +
    labs(title = "Reactome enrichment – Non-inflamed-specific genes")
  
  print(p_enrichment_noninflamed)
  
  ggsave(
    file.path(output_figures, "enrichment_noninflamed.png"),
    p_enrichment_noninflamed,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  message("No significant Reactome pathway found for non-inflamed-specific genes.")
}


# ============================================================
# 17. FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("NETWORK ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat("\nGlobal network statistics:\n")
print(network_global_statistics)

cat("\nEdge comparison:\n")
print(edge_comparison_summary)

cat("\nNetwork similarity:\n")
print(network_similarity)

cat("\nResults saved in:\n")
cat(normalizePath(output_dir), "\n")

cat("\nFolder structure:\n")
cat("results/network_analysis/\n")
cat("  tables/     CSV tables (statistics, edges, enrichment)\n")
cat("  figures/    PNG figures (networks and enrichment plots)\n")
cat("  networks/   RDS graph and enrichment objects\n")

cat("\nReminder: condition-specific edges are exploratory differences,\n")
cat("not proof of statistically significant differential regulation.\n")
