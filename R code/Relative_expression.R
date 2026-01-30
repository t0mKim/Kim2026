############################################################
# Relative expression analysis (edgeR, TMM, log2-CPM)
# Public / Anonymized version
############################################################

library(edgeR)

## 1) Paths (relative only) ----
# Suggested structure:
# project/
# ├── data/
# │   └── expression/
# └── results/
#     └── tables/

infile  <- file.path("data", "expression", "gene_count_table.txt")
outfile <- file.path("results", "tables", "selected_genes_log2CPM_TMM.txt")

## 2) Load count table ----
expr_df <- read.table(
  infile,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Create unique gene identifier
expr_df$unique_gene <- paste(expr_df$gene_name, expr_df$gene_id, sep = "_")

## 3) Select count columns dynamically (no UT/ANP/SC) ----
# Assumes count columns start after annotation columns
# Adjust 'meta_cols' if your table structure changes
meta_cols <- c("gene_id", "gene_name", "unique_gene")
count_cols <- setdiff(colnames(expr_df), meta_cols)

count_matrix <- as.matrix(expr_df[, count_cols])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- expr_df$unique_gene

## 4) Define genes of interest (kept minimal, no explanation here) ----
genes_of_interest <- c(
  "Acp5",
  "Nfatc1",
  "Oscar",
  "Rufy4",
  "Ocstamp",
  "Ostf1"
)

matched_rows <- expr_df[expr_df$gene_name %in% genes_of_interest,
                        c("gene_name", "unique_gene")]

# Enforce gene order
matched_rows <- matched_rows[
  match(genes_of_interest, matched_rows$gene_name),
]
matched_rows <- matched_rows[!is.na(matched_rows$gene_name), ]

## 5) edgeR normalization ----
y <- DGEList(counts = count_matrix)
y <- calcNormFactors(y, method = "TMM")

log2_cpm <- cpm(y, log = TRUE, prior.count = 1)

final_cpm <- log2_cpm[matched_rows$unique_gene, , drop = FALSE]

# Safe row names
rownames(final_cpm) <- make.unique(matched_rows$gene_name)

## 6) Output ----
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

write.table(
  final_cpm,
  file = outfile,
  sep = "\t",
  quote = FALSE
)

print(final_cpm)