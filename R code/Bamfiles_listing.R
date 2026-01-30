############################################################
# RNA-seq Gene Expression Quantification & DESeq2 Analysis
# Public / Anonymized Version
############################################################

## 1. Package setup ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_pkgs <- c("Rsubread", "DESeq2", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(Rsubread)
library(DESeq2)
library(ggplot2)

## 2. Define project paths (relative paths only) ----
# Project structure assumption:
# project/
# ├── data/
# │   ├── bam/
# │   └── annotation/
# └── scripts/

bam_dir <- file.path("data", "bam")
annotation_file <- file.path("data", "annotation", "genome.gtf")

## 3. Load BAM files ----
bam_files <- list.files(
  path = bam_dir,
  pattern = "\\.bam$",
  full.names = TRUE
)

stopifnot(length(bam_files) > 0)
print(bam_files)

## 4. Gene-level read counting using featureCounts ----
fc <- featureCounts(
  files = bam_files,
  annot.ext = annotation_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE
)

counts <- fc$counts
head(counts)

## 5. Sample metadata (anonymized) ----
# Replace with your own metadata file if needed
sample_names <- colnames(counts)

conditions <- factor(
  rep(c("Control", "Group1", "Group2"), each = 3),
  levels = c("Control", "Group1", "Group2")
)

colData <- data.frame(
  condition = conditions,
  row.names = sample_names
)

print(colData)

## 6. DESeq2 analysis ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# Basic filtering
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)
head(assay(vsd))

## 7. Session info for reproducibility ----
writeLines(
  capture.output(sessionInfo()),
  con = "sessionInfo.txt"
)