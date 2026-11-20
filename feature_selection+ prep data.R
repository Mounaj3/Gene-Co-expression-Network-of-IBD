
#========================================================================#
#  Saving the final datasets in csv format
#  Input : final_exp : expression matrix after preprocessing 
#          meta : data.table of metadata 
#          literature_gene_ids : vector of literature genes and Ensembl IDs
#          raw_counts : data.table 
#
#  Output: results/ : directory containing the final csv files
#          final_expression_top_MI_plus_literature :  final expression matrix
#          final_input_data.csv : expression data + metadata for 
#          final_selected_genes.rds : list of selected genes
#========================================================================#

library(data.table)
library(infotheo)
library(edgeR)

# ------------------------------------------------------------
# 1. Load processed expression and metadata
# ------------------------------------------------------------

n_top= 40
expr = fread("data/processed/expr_data.csv")
meta = fread("data/processed/meta_selected_samples.csv")


X = as.data.frame(expr[, -"samples_name"])
rownames(X) = expr$samples_name
X[] = lapply(X, as.numeric)

cat("Initial expression matrix:", nrow(X), "samples x", ncol(X), "genes\n")

# ------------------------------------------------------------
# Calculate mutual information matrix
# ------------------------------------------------------------

#Discretization
X_disc = discretize(as.data.frame(X), disc = "equalfreq",nbins = 10 )

# compute MI 
# MI takes around 3H30 to run for 3000 genes
# MI matrix is available 
# MI=fread("results/MI_matrix.csv")
# MI=as.matrix(MI[,-1], rownames = MI$V1)

t0 = Sys.time()
#MI = mutinformation(X_disc)
diag(MI) = 0

cat("MI calculation completed in",round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),"minutes\n")

rownames(MI) = colnames(X_disc)
colnames(MI) = colnames(X_disc)


# Save MI matrix
fwrite(
  as.data.table(MI, keep.rownames = "gene"),
  "results/MI_matrix.csv"
)


# ------------------------------------------------------------
#  Rank genes
# ------------------------------------------------------------

# Mean MI with all other genes
mean_MI = rowSums(MI) / (ncol(MI) - 1) 
# Decreasing order   
mean_MI = sort(mean_MI, decreasing = TRUE)

MI_ranking = data.frame(
  gene = names(mean_MI),
  mean_MI = as.numeric(mean_MI)
)

# Save complete ranking
fwrite(
  MI_ranking,
  "results/MI_gene_ranking.csv"
)

# Select top 40 genes
top_MI = names(mean_MI)[
  seq_len(min(n_top, length(mean_MI)))
]

cat("Top 40 genes selected by mean MI\n", top_MI)

# ------------------------------------------------------------
# Literature-based genes 
# ------------------------------------------------------------

literature_gene_ids  = c(
  GBP1 = "ENSG00000117228",
  TAP1 = "ENSG00000168394",
  `HLA-E` = "ENSG00000204592",
  `HLA-B` = "ENSG00000234745",
  STAT1 = "ENSG00000115415",
  PARP9 = "ENSG00000138496",
  PARP14 = "ENSG00000173193",
  NFKB1 = "ENSG00000109320",
  STAT3 = "ENSG00000168610",
  JUN = "ENSG00000177606",
  FOS = "ENSG00000170345",
  SMAD2 = "ENSG00000175387"
)

# ------------------------------------------------------------
# Match literature genes to expression data
# ------------------------------------------------------------


expression_gene_ids = colnames(X)

# Identify literature genes already present in X
lit_positions = match(literature_gene_ids, expression_gene_ids)
present = !is.na(lit_positions)

cat(
  "Literature genes already present:",
  sum(present), "/", length(literature_gene_ids), "\n"
)

cat("Literature genes from X (commun):\n")
print(names(literature_gene_ids)[present])
cat("Literature genes missing from X:\n")
print(names(literature_gene_ids)[!present])


# ------------------------------------------------------------
# Retrieve missing genes from the original raw count data
# ------------------------------------------------------------
if (any(!present)) {
  
  # Read original raw count matrix
  raw_counts = fread("data/row/GSE193677_MSCCR_Biopsy_counts.txt")
  
  # Gene IDs (first column)
  raw_gene_ids = raw_counts$V1
  
  # Combine top MI genes + literature genes
  final_genes = unique(c(top_MI, literature_gene_ids))

  
  # Identify all selected genes missing from X
  expression_gene_ids = colnames(X)
  missing_genes = setdiff(final_genes, expression_gene_ids)
  
  # Find missing genes in raw counts
  raw_positions = match(missing_genes, raw_gene_ids)
  
  if (anyNA(raw_positions)) {
    stop("Some selected genes were not found in the raw count file")
  }
  
  # Extract raw counts
  raw_count_matrix = as.matrix(
    raw_counts[raw_positions, -1, with = FALSE]
  )
  
  rownames(raw_count_matrix) = missing_genes
  storage.mode(raw_count_matrix) = "numeric"
  
  # Keep only selected samples, in the same order as X
  raw_count_matrix = raw_count_matrix[
    , rownames(X), drop = FALSE
  ]
  
  # TMM normalization + log2-CPM
  y = DGEList(counts = raw_count_matrix)
  y = calcNormFactors(y)
  logCPM_full = cpm(y, log = TRUE)
  
  # Convert to samples x genes
  recovered_genes = t(logCPM_full)
  
  # Add recovered genes to X
  X = cbind(X, recovered_genes)
}

# ------------------------------------------------------------
# Combine top 40 MI genes + all literature genes
# ------------------------------------------------------------

expression_gene_ids =  colnames(X)

final_genes = unique(c(top_MI, literature_gene_ids))


# Check whether all selected genes are present
missing_genes = setdiff(final_genes, expression_gene_ids)

print(missing_genes)

if (length(missing_genes) > 0) {
  stop(
    "Some selected genes are still absent: ",
    paste(missing_genes, collapse = ", ")
  )
}

# Final expression matrix
final_gene_names = colnames(X)[
  match(final_genes, expression_gene_ids)
]

final_exp = X[, final_gene_names, drop = FALSE]

cat("Final number of selected genes:", ncol(final_exp), "\n")

# Final expression expr
final_gene_names = colnames(X)[match(final_genes, expression_gene_ids)]
final_exp = X[, final_gene_names,drop = FALSE]

cat("Top MI genes:", length(unique(top_MI)), "\n")
cat("Literature genes:", length(literature_gene_ids), "\n")
cat("Final unique genes:", ncol(final_exp), "\n")


# ------------------------------------------------------------
# Build final data
# ------------------------------------------------------------

# Align metadata to expression sample order
meta = as.data.frame(meta)
meta = meta[match(rownames(X), meta$samples_name),,drop = FALSE]
stopifnot(identical(meta$samples_name, rownames(X)))

meta_variables = data.frame(
  Disease = meta$ibd_disease,                # CD / UC
  Region  = meta$region_groups,              # Ileum / Colon / Rectum
  Sex     = meta$demographics_gender,        # Female / Male
  EndoSev = meta$ibd_endoseverity_4levels,   # Inactive < Mild < Moderate < Severe
  CRP     = meta$crp_jjmgl_log2,             # systemic inflammation (log2)
  Histo   = meta$ghas_sum7                   # histology of the biopsy (many NA)
)


# Expression: samples x selected genes
final_data = data.frame(
  samples_name = rownames(final_exp),
  final_exp,
  meta_variables,
  check.names = FALSE
)
dim(final_data)
cat(
  "Final MIIC dataset:",
  nrow(final_data), "samples x",
  ncol(final_data) - 1, "variables + sample ID\n"
)



# Save final datasets
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Final expression matrix
fwrite( data.frame(samples_name = rownames(final_exp),final_exp),
  "results/final_expression_top_MI_plus_literature.csv"
)

# Final MIIC dataset: expression + metadata
fwrite(final_data,"results/final_input_data.csv")

# Save selected gene list
fwrite(data.frame(gene = final_genes),"results/final_selected_genes.csv")

# Save literature genes that were not found
fwrite(
  data.frame(
    gene_symbol = names(lit_ens)[!present],
    Ensembl_ID = unname(lit_ens[!present])
  ),
  "results/missing_literature_genes.csv"
)


cat("===============================\nFinal results\n")
cat("Top MI genes:", length(top_MI), "\n")
cat("Literature genes found:", sum(present), "\n")
cat("Final selected genes:", length(final_genes), "\n")
cat("Samples:", nrow(final_data), "\n")

