
############################################################
# Exploratory Data Analysis (EDA)
# Overview of gene expression variability, sample correlations,
# expression patterns, PCA and mutual information between genes.
############################################################

library(data.table)
library(BioNERO)
library(matrixStats)
library(ggplot2)
library(ComplexHeatmap)

# 1. Load data
expr = fread("data/processed/expr_data.csv")
meta = fread("data/processed/meta_selected_samples.csv")

plot_dir = "results/EDA"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Prepare expression matrix (genes x samples)
x = as.matrix(t(as.matrix(expr[, -"samples_name"])))
colnames(x) = expr$samples_name
storage.mode(x) = "numeric"

col_meta = as.data.frame(meta)
rownames(col_meta) = col_meta$samples_name

if (any(!colnames(x) %in% rownames(col_meta))) {
  stop("Some samples are missing from the metadata.")
}

col_meta = col_meta[
  colnames(x),
  c("typere", "ibd_disease", "region_groups"),
  drop = FALSE
]
stopifnot(identical(colnames(x), rownames(col_meta)))
col_meta[] = lapply(col_meta, as.factor)


# 4. Gene variability
variance_values = rowVars(x, na.rm = TRUE)
variance_values = sort(
  variance_values[is.finite(variance_values)],
  decreasing = TRUE
)

variance_df = data.frame(
  gene = names(variance_values),
  Variance = as.numeric(variance_values)
)

# Plot 1: Distribution of gene expression variance
p1 = ggplot(variance_df, aes(x = Variance)) +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(
    title = "Distribution of gene expression variance",
    x = "Variance",
    y = "Number of genes"
  )


ggsave(
  file.path(plot_dir, "01_variance_distribution.png"),
  p1, width = 7, height = 5
)

# Plot 2: Top 30 most variable genes
top30 = head(variance_df, 30)

p2 = ggplot(top30, aes(x = reorder(gene, Variance), y = Variance)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top 30 most variable genes",
    x = "Gene",
    y = "Variance"
  )

ggsave(
  file.path(plot_dir, "02_top30_variable_genes.png"),
  p2, width = 7, height = 8
)

# Plot 3: Sample-to-sample correlation using ALL genes
p3 = plot_heatmap(
  x,
  col_metadata = col_meta,
  type = "samplecor",
  cor_method = "pearson",
  show_rownames = FALSE,
  show_colnames = FALSE
)

png(
  filename = file.path(plot_dir, "03_sample_correlation_all_genes.png"),
  width = 2000,
  height = 1800,
  res = 300
)

ComplexHeatmap::draw(p3)
dev.off()

# Plot 4: Expression heatmap using ALL genes
p4 = plot_heatmap(
  x,
  col_metadata = col_meta,
  type = "expr",
  show_rownames = FALSE,
  show_colnames = FALSE
)

png(
  filename = file.path(plot_dir, "04_expression_heatmap_genes.png"),
  width = 2400,
  height = 2000,
  res = 300
)

ComplexHeatmap::draw(p4)
dev.off()

# Plot 5: Expression heatmap using the TOP 50 most variable genes
top50_genes = head(variance_df$gene, min(50, nrow(variance_df)))

x_top50 = x[top50_genes, , drop = FALSE]

p5 = plot_heatmap(
  x_top50,
  col_metadata = col_meta,
  type = "expr",
  show_rownames = TRUE,
  show_colnames = FALSE
)

png(
  filename = file.path(plot_dir, "05_expression_heatmap_top50_genes.png"),
  width = 2400,
  height = 2000,
  res = 300
)

ComplexHeatmap::draw(p5)
dev.off()

# Plot 6: PCA according to disease status and region
p6= plot_PCA(
  x,
  metadata = col_meta[, c("region_groups", "ibd_disease"), drop = FALSE]
)

ggsave(
  file.path(plot_dir, "06_PCA_disease_region.png"),
  p6, width = 8, height = 6
)


# Plot 7: Mutual Information (MI) heatmap

MI = fread("results/MI_matrix.csv")
MI = as.matrix(MI[, -1, with = FALSE])
rownames(MI) = fread("results/MI_matrix.csv")$V1
storage.mode(MI) = "numeric"
MI[1:3,1:3]
# Save MI heatmap
png(
  filename = file.path(plot_dir, "07_MI_heatmap.png"),
  width = 2400,
  height = 2200,
  res = 300
)

ComplexHeatmap::draw(
  ComplexHeatmap::Heatmap(
    MI[1:3000,1:3000],
    name = "Mutual Information",
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = "Mutual Information between genes",
    row_title = "Genes"

))


dev.off()

# Plot 8: Expression heatmap of the 53 selected genes


final_x = fread("results/final_expression_top_MI_plus_literature.csv")

x_selected = as.matrix(t(as.matrix(final_x[, -"samples_name"])))
colnames(x_selected) = final_x$samples_name
storage.mode(x_selected) = "numeric"
selected_meta = col_meta[colnames(x_selected), , drop = FALSE]


# Generate and save heatmap
p_selected = plot_heatmap(
  x_selected,
  col_metadata = selected_meta,
  type = "expr",
  show_rownames = FALSE,
  show_colnames = FALSE
)


png(
  filename = file.path(plot_dir, "08_expression_heatmap_53_selected_genes.png"),
  width = 2400,
  height = 2000,
  res = 300
)

ComplexHeatmap::draw(p_selected)
dev.off()
