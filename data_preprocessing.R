#========================================================================#
#  Preprocessing of the MSCCR expression data (GSE193677)
#  Input : expr = data.table ou RAW COUNTS
#                 samples in rows 
#                 "samples_name"  column  + one column per gene (Ensembl IDXXX)
#  Output: expr_clean = same format, log2-CPM (TMM), filtered genes
#===========================================================================#

library(data.table)
library(dplyr)
library(BioNERO)
library(edgeR)
library(matrixStats)

set.seed(123)

#---------------------------------------------------------------
# Keep known genes (protein-coding,long non-coding RNAs : lncRNA)
#
#---------------------------------------------------------------
filter_gene_type = function(dt, annot, keep_types = c("protein-coding")) {
  annot = annot[EnsemblGeneID != ""]
  annot = annot[!duplicated(EnsemblGeneID)]
  
  genes = setdiff(names(dt), "samples_name")          
  info  = annot[match(genes, EnsemblGeneID)]
  
  cat("Gene types in the matrix:\n")
  print(table(info$GeneType, useNA = "ifany"))
  
  keep = genes[!is.na(info$GeneType) & info$GeneType %in% keep_types]
  cat("Kept:", length(keep), "genes of type", paste(keep_types, collapse = ", "), "\n")
  dt[, c("samples_name", keep), with = FALSE]
}

#---------------------------------------------------------------
# Remove mitochondrial and ribosomal genes
# NC_012920.1 = RefSeq accession of the human mitochondrial genome
# Ribosom
#---------------------------------------------------------------
remove_ribo_mito = function(dt, annot) {
  annot = annot[EnsemblGeneID != ""][!duplicated(EnsemblGeneID)]
  genes = setdiff(names(dt), "samples_name")
  info  = annot[match(genes, EnsemblGeneID)]       
  
  is_mito = info$ChrAcc == "NC_012920.1"
  is_ribo = (grepl("^RP[LS][0-9]|^RPLP[0-9]|^MRP[LS][0-9]", info$Symbol) |
                info$Symbol %in% c("FAU", "UBA52")) &
    !grepl("^RPS6K", info$Symbol)     # RPS6K = kinases, keep them
  is_mito[is.na(is_mito)] = FALSE
  is_ribo[is.na(is_ribo)] = FALSE
  
  cat("Mitochondrial genes removed:", sum(is_mito), "\n")
  cat("Ribosomal genes removed:    ", sum(is_ribo & !is_mito), "\n")
  
  keep = genes[!(is_mito | is_ribo)]
  cat("Genes left:                 ", length(keep), "\n")
  dt[, c("samples_name", keep), with = FALSE]
}

#---------------------------------------------------------------
#  Missing values, non-expressed genes, normalisation,
#    constant genes, optional top_n and outliers (BioNERO + edgeR)
#---------------------------------------------------------------
preprocess_bionero = function(dt,
                               min_exp = 20,           # raw-count threshold (reads)
                               min_prop = 0.4,         # fraction of samples >= min_exp
                               top_n = NULL,           # n most variable genes (NULL = all)
                               remove_outliers = FALSE) {
  
  # BioNERO format (genes in rows, samples in columns)
  exp = as.data.frame(t(as.matrix(dt[, -"samples_name"])))
  colnames(exp) = dt$samples_name
  cat("Start:", nrow(exp), "genes,", ncol(exp), "samples\n")
  
  # Missing values -> 0
  exp = replace_na(exp)
  
  #  Non-expressed genes 
  exp = remove_nonexp(exp, method = "percentage",
                       min_exp = min_exp, min_percentage_samples = min_prop)
  cat("After non-expressed filter:", nrow(exp), "genes\n")
  
  # Normalisation: TMM + log2-CPM (edgeR)
  y   = calcNormFactors(DGEList(counts = as.matrix(exp)))
  exp = as.data.frame(cpm(y, log = TRUE))
  
  # Constant genes (zero variance)
  v = rowVars(as.matrix(exp))
  cat("Constant genes removed:", sum(v == 0), "\n")
  exp = exp[v > 0, ]
  
  # Optional: keep the top_n most variable genes
  if (!is.null(top_n)) {
    exp = filter_by_variance(exp, n = top_n)
    cat("After variance filter:", nrow(exp), "genes\n")
  }
  
  # Optional: outlier samples (Z.K method)
  if (remove_outliers) {
    exp = ZKfiltering(exp, cor_method = "pearson")
  }
  
  cat("End:", nrow(exp), "genes,", ncol(exp), "samples\n")
  

  data.table(samples_name = colnames(exp), t(as.matrix(exp)))
}

#==============================================================#
#  Run the pipeline (on RAW counts)
#==============================================================#

genes_annotation=fread("data/selected_genes_annotation.csv")
expr=fread("data/counts_selected_samples.csv")

# 1. Known genes: protein-coding + long non-coding RNAs
expr_selected_genes = filter_gene_type(expr, genes_annotation,
                                        keep_types = c("protein-coding", "ncRNA"))

# 2. Remove mitochondrial and ribosomal genes
expr_selected_genes = remove_ribo_mito(expr_selected_genes, genes_annotation)

# 3. Missing values, non-expressed, normalisation, constant genes
expr_clean = preprocess_bionero(expr_selected_genes, min_exp = 20, min_prop = 0.4)

# Annotation of the genes kept (in the same order of counts matrix)
annotation_kept = genes_annotation[match(names(expr_clean)[-1], EnsemblGeneID)]

# Save
dir.create("data", showWarnings = FALSE)
fwrite(expr_clean, "data/expr_clean_logcpm.csv")
fwrite(annotation_kept, "data/selected_genes_annotation.csv")

p = plot_heatmap(expr_clean, type = "samplecor", show_rownames = FALSE)
