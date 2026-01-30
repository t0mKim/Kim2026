############################################################
# Venn diagram
# Public / Anonymized version
############################################################

# Ensure counts is a matrix
counts_mat <- as.matrix(counts)

# Use anonymized condition labels
# Assumes you have: colData with rownames matching colnames(counts_mat)
stopifnot(all(colnames(counts_mat) %in% rownames(colData)))

# Align metadata order to count matrix columns
meta <- colData[colnames(counts_mat), , drop = FALSE]
cond <- as.character(meta$condition)

# Define two groups to compare (anonymized)
group_A_label <- "Control"   # e.g., Control
group_B_label <- "Group2"    # e.g., Treatment-like group

group_A_samples <- colnames(counts_mat)[cond == group_A_label]
group_B_samples <- colnames(counts_mat)[cond == group_B_label]

stopifnot(length(group_A_samples) > 0, length(group_B_samples) > 0)

# Expression thresholding rule:
# gene is "expressed" if counts >= threshold in at least half of samples in each group
threshold <- 10
min_n_A <- ceiling(length(group_A_samples) / 2)
min_n_B <- ceiling(length(group_B_samples) / 2)

group_A_expr_genes <- rownames(counts_mat)[
  rowSums(counts_mat[, group_A_samples, drop = FALSE] >= threshold) >= min_n_A
]
group_B_expr_genes <- rownames(counts_mat)[
  rowSums(counts_mat[, group_B_samples, drop = FALSE] >= threshold) >= min_n_B
]

# Packages
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(grid)
library(VennDiagram)

# Create the Venn diagram (generic labels, generic title)
venn_plot <- venn.diagram(
  x = list(
    GroupA = group_A_expr_genes,
    GroupB = group_B_expr_genes
  ),
  filename = NULL,
  fill = c("grey70", "grey40"),
  alpha = c(0.5, 0.5),
  cat.cex = 1.3,
  cex = 1.2,
  main = "Venn Diagram: Expressed Genes (GroupA vs GroupB)",
  scaled = FALSE
)

grid.draw(venn_plot)