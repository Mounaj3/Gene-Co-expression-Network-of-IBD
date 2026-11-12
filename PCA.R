library(data.table)
library(dplyr)
library(ggplot2)
library(corrplot)
library(matrixStats)

#==============================================================#
#  Feature Selection:find features that drives most the variance
#  PCA - exploratory analysis
#  Input : expr = data.table 
#          samples in rows, 
#         "samples_name+ genes, log2-CPM) ; meta = sample metadata
#          
#==============================================================#

#---------------------------------------------------------------
# Matrix samples x genes + metadata 
#---------------------------------------------------------------
expr=fread( "data/expr_clean_logcpm.csv")
meta =fread("data/meta_selected_samples.csv")

#setdiff(meta$samples_name, expr$samples_name)
str(expr)
str(meta)
samples_label=meta[,c("samples_name","typere","region_groups","ibd_disease")]
samples_label=samples_label[match(expr$samples_name, meta$samples_name), ]

#df=merge(as.data.table(samples_label), as.data.table(expr), by="samples_name")

X = as.matrix(expr[, -"samples_name"])
rownames(X) = expr$samples_name

meta_pca = as.data.frame(meta)[match(rownames(X), meta$samples_name), ]
stopifnot(all(meta_pca$samples_name == rownames(X))) 


#---------------------------------------------------------------
# Run PCA on the 2,000 most variable genes ( impossible on all)
#---------------------------------------------------------------
top = order(colVars(X), decreasing = TRUE)[1:min(2000, ncol(X))]
pca = prcomp(X[, top], center = TRUE, scale. = TRUE)



# sum of squared distance between genes and pc (eigenvalues)
#--> allow as to determine proportion of variation of pc
# pc 1 accounts for 7 
pca$rotation[1:3,1:3]
head(pca$sdev)
head(sort(abs(pca$rotation[, 1]), decreasing = TRUE), 20)
var_expl = 100 * pca$sdev^2 / sum(pca$sdev^2) # % variance per PC



#---------------------------------------------------------------
#  Plot: variance explained by the first 10 PCs
#---------------------------------------------------------------
ggplot(data.frame(PC = 1:10, var = var_expl[1:10]), aes(PC, var)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(var, 1)), vjust = -0.3, size = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Variance explained by each PC", x = "Principal component",
       y = "Variance explained (%)") +
  theme_minimal()
#  Top 20 genes driving PC1 
plot(var_expl[1:20], type = "b", pch = 19,
     xlab = "PC", ylab = "Proportion of variance", main = "Scree plot")

#---------------------------------------------------------------
#  Gene selection 
#---------------------------------------------------------------

#  set K = the number of PCs before the "elbow"
K = 5 # Most representative PCs of our data 
# ---- Get the loadings (genes x PCs) ----

loadings = abs(pca$rotation[, 1:K])  

# ----  Take the top genes of each PC ----
n_top = 150 # genes kept per PC

top_genes = c()
for (k in 1:K) {
  genes_k   = names(sort(loadings[, k], decreasing = TRUE))[1:n_top]
  top_genes = c(top_genes, genes_k)
}
selected_genes = unique(top_genes)     # remove duplicate



#---------------------------------------------------------------
#  PC1 vs PC2 coloured by each tissus status, region and IBD type
#---------------------------------------------------------------
scores = cbind(meta_pca, as.data.frame(pca$x[, 1:5]))

plot_pc = function(var, title) {
  ggplot(scores, aes(PC1, PC2, colour = .data[[var]])) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_brewer(palette = "Set2") +
    labs(title = title,
         x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
         y = paste0("PC2 (", round(var_expl[2], 1), "%)"),
         colour = "") +
    theme_minimal()
}
plot_pc("region_groups",     "PCA coloured by region")
plot_pc("typere",      "PCA coloured by inflammation status")
plot_pc("ibd_disease", "PCA coloured by diagnosis")

