############################################################
# Running Enrichment Score (pre-ranked; GO-based geneset)
# Public / Anonymized version
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

## ---- Inputs (relative paths only) ----
ranked_file <- file.path("data", "ranked_gene_list.tsv")
out_dir     <- file.path("results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Keep the identifier configurable; avoid semantic labels in code/comments
geneset_id <- "GO:0050728"

# Optional: anonymized tag for output filenames (so files don't reveal topic)
geneset_tag <- "GENESET_A"

## ---- 1) Load ranked list ----
df <- read_tsv(ranked_file, show_col_types = FALSE) |> as.data.frame()

stopifnot(ncol(df) >= 3)

# Standardize expected columns without exposing original meanings
colnames(df)[1] <- "gene_full"   # e.g., "Symbol(EnsemblID...)"
colnames(df)[3] <- "metric"

# Extract SYMBOL from "Gene(....)" format
df$gene <- sub("\\(.*", "", df$gene_full)

# Keep needed columns and clean
df <- df[, c("gene", "metric")]
df$metric <- as.numeric(df$metric)
df <- df[!is.na(df$metric) & !is.na(df$gene) & df$gene != "", , drop = FALSE]

# Handle duplicate symbols (one common rule: take max score)
df <- aggregate(metric ~ gene, data = df, FUN = max)

# Sort descending (pre-ranked)
df <- df[order(-df$metric), ]
rownames(df) <- NULL

stats <- setNames(df$metric, df$gene)
N <- length(stats)

## ---- 2) Fetch geneset (SYMBOL) from GO ID ----
sym_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = geneset_id,
  keytype = "GO",
  columns = c("SYMBOL", "ENTREZID")
)

geneset <- unique(na.omit(sym_map$SYMBOL))

# Overlap check
common <- intersect(names(stats), geneset)
message("Overlap: ", length(common), " of ", length(geneset))

tag <- names(stats) %in% geneset
if (sum(tag) == 0) stop("No overlap between ranked list and selected geneset.")

normHit  <- sum(abs(stats[tag]))   # weighted by |metric|, p=1
normMiss <- N - sum(tag)

if (normHit == 0)  stop("normHit is 0 (all overlapping metrics are zero).")
if (normMiss == 0) stop("normMiss is 0 (geneset covers all ranked genes).")

## ---- 3) Running ES (weighted) ----
runningES <- numeric(N)
cum <- 0
for (i in seq_len(N)) {
  cum <- cum + if (tag[i]) abs(stats[i]) / normHit else -1 / normMiss
  runningES[i] <- cum
}

## ---- 4) Outputs ----
output_all <- data.frame(
  gene      = names(stats),
  rank      = seq_len(N),
  runningES = runningES,
  metric    = as.numeric(stats),
  is_hit    = tag,
  stringsAsFactors = FALSE
)

write.table(
  output_all,
  file.path(out_dir, paste0("RUNNING_ES_", geneset_tag, "_ALL.tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

output_hits <- output_all[output_all$is_hit, , drop = FALSE]
write.table(
  output_hits,
  file.path(out_dir, paste0("RUNNING_ES_", geneset_tag, "_HITS.tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

print(head(output_all, 3))
print(tail(output_all, 3))