############################################################
# Running Enrichment Score (pre-ranked; GO-based geneset)
# Public / Anonymized version
############################################################

# ---------- Packages ----------
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ---------- I/O (GitHub-safe, low-context) ----------
infile  <- file.path("data", "enrichment", "term_results.xlsx")
outfile <- file.path("results", "figures", "dotplot_terms.png")
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

# ---------- Load data ----------
df <- read_excel(infile)

# ---------- GeneRatio to numeric ----------
if (is.character(df$GeneRatio)) {
  gr <- gsub("\\s+", "", df$GeneRatio)  # remove spaces
  df$GeneRatioNum <- sapply(strsplit(gr, "/"), function(x) {
    if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else as.numeric(x[1])
  })
} else {
  df$GeneRatioNum <- as.numeric(df$GeneRatio)
}

# ---------- padj numeric + clean ----------
df$padj <- suppressWarnings(as.numeric(df$padj))
df <- subset(df, !is.na(padj) & !is.na(GeneRatioNum) & !is.na(Count) & !is.na(Description))

# (Optional) keep only significant terms:
# df <- subset(df, padj <= 0.05)

# Order y by padj (small on top)
df$Description <- fct_rev(fct_reorder(df$Description, df$padj))

# ---------- Plot ----------
p <- ggplot(df, aes(x = GeneRatioNum, y = Description, size = Count, color = padj)) +
  geom_point() +
  
  # Color scale (cap at 0.05 if desired)
  scale_color_gradient(
    low    = "red",
    high   = "blue",
    limits = c(0, 0.05),
    oob    = squish,
    breaks = c(0, 0.025, 0.05),
    labels = c("0", "0.025", "0.05")
  ) +
  
  # Size scale (automatic breaks)
  scale_size_continuous(range = c(3, 8)) +
  
  # Zoom without dropping data
  coord_cartesian(xlim = c(0, 0.10)) +
  
  labs(
    x = "Gene Ratio",
    y = "Pathway/Term",
    color = "Adjusted p-value",
    size  = "Gene Count"
  ) +
  theme_bw() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    axis.text.y       = element_text(size = 10)
  )

ggsave(outfile, plot = p, width = 7.5, height = 5.0, dpi = 300, units = "in", bg = "transparent")