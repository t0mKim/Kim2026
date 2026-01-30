############################################################
# ComplexHeatmap 
# Public / Anonymized version
############################################################

# Packages
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) install.packages("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

library(ComplexHeatmap)
library(circlize)

## 1) Paths (relative only) ----
# Suggested structure:
# project/
# ├── data/
# │   └── heatmap/
# └── results/
#     └── figures/

infile  <- file.path("data", "heatmap", "heatmap_input.csv")
outfile <- file.path("results", "figures", "heatmap.png")

## 2) Load matrix ----
heatmap_mat <- read.csv(
  infile,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
heatmap_mat <- as.matrix(heatmap_mat)

## 3) Optional: drop a group by pattern (anonymized) ----
# If you need to drop a subset of samples, use generic patterns like "Group1"
# heatmap_mat <- heatmap_mat[, !grepl("^Group1", colnames(heatmap_mat)), drop = FALSE]

## 4) QC filtering ----
heatmap_mat <- heatmap_mat[rowSums(heatmap_mat) > 0, , drop = FALSE]          # remove all-zero genes
heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, var) > 0, , drop = FALSE]    # remove zero-variance genes

## 5) Row-wise z-score ----
scaled_mat <- t(scale(t(heatmap_mat)))

## 6) Anonymize sample names (no UT/ANP etc.) ----
# Keep column order but rename generically as Sample01, Sample02, ...
colnames(scaled_mat) <- sprintf("Sample%02d", seq_len(ncol(scaled_mat)))

## 7) Clustering ----
hc <- hclust(dist(scaled_mat, method = "euclidean"), method = "ward.D2")
dend <- as.dendrogram(hc)

## 8) Color mapping (generic) ----
# NOTE: Avoid condition-signaling color schemes if you want extra anonymity.
# Keeping a standard diverging mapping around 0 is fine.
color_mapping <- colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("green","darkgreen","black","darkred","red")
)

ht <- Heatmap(
  scaled_mat,
  name = "Z-score",
  cluster_rows = dend,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = color_mapping
)

## 9) Save figure (relative path, no private info) ----
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

png(outfile, bg = "transparent", width = 1300, height = 2500, res = 300)
draw(ht, background = "transparent")
dev.off()