
# =============================================================================
# PC INFERENCE
# =============================================================================

library(data.table)
library(pcalg)
library(igraph)

# 1. PARAMETRES ---------------------------------------------------------------

input_file <- "results/results/final_input_data.csv"
output_dir <- "results/PC"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.01
m.max <- 2

metadata_cols <- c(
  "typere", "Disease", "Region", "Sex", "EndoSev", "CRP"
)

data <- fread(input_file)
gene_cols <- setdiff(names(data), metadata_cols)


# 2. FONCTION D'INFERENCE PC ---------------------------------------------------

run_PC <- function(X) {
  
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  if (anyNA(X) || any(!is.finite(X)) ||
      any(apply(X, 2, sd) == 0)) {
    stop("Invalid data: missing, non-finite or constant values.")
  }
  
  pcalg::pc(
    suffStat = list(C = cor(X), n = nrow(X)),
    indepTest = gaussCItest,
    alpha = alpha,
    labels = colnames(X),
    skel.method = "stable",
    m.max = m.max,
    solve.confl = TRUE
  )
}


# 3. EXTRACTION DES ARETES ----------------------------------------------------

extract_edges <- function(pc_fit) {
  
  amat <- as(pc_fit, "amat")
  genes <- colnames(amat)
  
  ij <- which(
    upper.tri(amat) & (amat != 0 | t(amat) != 0),
    arr.ind = TRUE
  )
  
  if (nrow(ij) == 0) {
    return(data.table(
      from = character(),
      to = character(),
      type = character()
    ))
  }
  
  i <- ij[, 1]
  j <- ij[, 2]
  
  edges <- data.table(
    from = genes[i],
    to = genes[j],
    type = ifelse(
      amat[cbind(i, j)] == amat[cbind(j, i)],
      "undirected",
      "directed"
    )
  )
  
  reverse <- amat[cbind(i, j)] == 1 &
    amat[cbind(j, i)] == 0
  
  edges[reverse, `:=`(from = to, to = from)]
  
  edges
}


# 4. DISPERSION DES NOEUDS ----------------------------------------------------

repel_layout <- function(coords, min_dist = 0.16, iterations = 500) {
  
  n <- nrow(coords)
  if (n <= 1) return(coords)
  
  coords <- scale(coords, center = TRUE, scale = FALSE)
  coords <- coords / max(abs(coords))
  
  for (k in seq_len(iterations)) {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        
        delta <- coords[j, ] - coords[i, ]
        distance <- sqrt(sum(delta^2))
        
        if (distance < 1e-10) {
          delta <- runif(2, -1e-5, 1e-5)
          distance <- sqrt(sum(delta^2))
        }
        
        if (distance < min_dist) {
          displacement <- (min_dist - distance) / 2
          coords[i, ] <- coords[i, ] - displacement * delta / distance
          coords[j, ] <- coords[j, ] + displacement * delta / distance
        }
      }
    }
  }
  
  scale(coords, center = TRUE, scale = FALSE)
}


# 5. EXPORT CSV ET PNG --------------------------------------------------------

export_network <- function(edges, genes, name, node_color) {
  
  # CSV des aretes
  fwrite(
    edges,
    file.path(output_dir, paste0("PC_", name, "_edges.csv"))
  )
  
  # Graphe
  graph <- graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = genes)
  )
  
  # Taille des noeuds proportionnelle au degre
  degree <- rowSums(
    (amat != 0) | (t(amat) != 0)
  )
  
  node_sizes <- if (max(degree) == min(degree)) {
    rep(48, length(degree))
  } else {
    36 + 24 * (degree - min(degree)) / (max(degree) - min(degree))
  }
  
  V(graph)$size <- node_sizes
  
  # Style des aretes
  E(graph)$color <- ifelse(
    E(graph)$type == "directed", "#263746", "#8795A1"
  )
  E(graph)$width <- ifelse(
    E(graph)$type == "directed", 2.2, 1.6
  )
  E(graph)$arrow.mode <- ifelse(
    E(graph)$type == "directed", 2, 0
  )
  
  # Layout
  set.seed(123)
  
  coords <- layout_with_fr(graph, niter = 5000)
  coords <- repel_layout(coords) * 10
  
  padding <- pmax(
    apply(coords, 2, function(x) diff(range(x))) * 0.25,
    1
  )
  
  limits_x <- range(coords[, 1]) + c(-padding[1], padding[1])
  limits_y <- range(coords[, 2]) + c(-padding[2], padding[2])
  
  center <- colMeans(coords)
  
  angles <- atan2(
    coords[, 2] - center[2],
    coords[, 1] - center[1]
  )
  
  # PNG
  png(
    file.path(output_dir, paste0("PC_network_", name, ".png")),
    width = 4800, height = 4200, res = 300
  )
  
  par(mar = c(5, 3, 5, 3), bg = "white")
  
  plot(
    graph,
    layout = coords,
    rescale = FALSE,
    xlim = limits_x,
    ylim = limits_y,
    asp = 1,
    
    vertex.color = node_color,
    vertex.frame.color = "#34495E",
    vertex.frame.width = 1.1,
    vertex.label = V(graph)$name,
    vertex.label.color = "#17212B",
    vertex.label.cex = 0.65,
    vertex.label.font = 2,
    vertex.label.dist = 1.2,
    vertex.label.degree = angles,
    
    edge.color = E(graph)$color,
    edge.width = E(graph)$width,
    edge.arrow.mode = E(graph)$arrow.mode,
    edge.arrow.size = 0.65,
    edge.arrow.width = 1.1,
    
    vertex.size = node_sizes,
    main = paste("PC gene network -", name),
    cex.main = 1.5
  )
  
  legend(
    "bottom",
    legend = c(
      paste("Genes -", name),
      "Directed edge",
      "Undirected edge",
      "Node size = degree"
    ),
    pch = c(21, NA, NA, 21),
    pt.bg = c(node_color, NA, NA, "white"),
    col = c("#34495E", "#263746", "#8795A1", "#34495E"),
    lty = c(NA, 1, 1, NA),
    lwd = c(NA, 2.2, 1.6, NA),
    horiz = TRUE,
    bty = "n",
    cex = 0.9
  )
  
  dev.off()
}


# 6. EXECUTION DES TROIS RESEAUX ----------------------------------------------

run_network <- function(dataset, name, node_color) {
  
  cat("\n==========", toupper(name), "==========\n")
  
  X <- as.matrix(dataset[, ..gene_cols])
  storage.mode(X) <- "double"
  
  if (nrow(X) < 10) stop("Not enough samples: ", name)
  
  pc_fit <- run_PC(X)
  edges <- extract_edges(pc_fit)
  
  # Matrice d'adjacence necessaire au calcul des degres
  amat <- as(pc_fit, "amat")
  
  export_network(edges, colnames(X), name, node_color)
  
  cat("Samples:", nrow(X), "\n")
  cat("Genes:", ncol(X), "\n")
  cat("Edges:", nrow(edges), "\n")
}

# Global
run_network(
  data,
  "global",
  "#83B5D8"
)

# Inflammé
run_network(
  data[typere == "I"],
  "inflamed",
  "#E99A91"
)

# Non inflammé
run_network(
  data[typere == "NonI"],
  "noninflamed",
  "#83BDA5"
)

cat("\nPC ANALYSIS COMPLETED\nResults saved in:", output_dir, "\n")
