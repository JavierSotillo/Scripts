# =============================================================================
# Volcano plot from a table of fold changes and adjusted p-values (Q-values)
#
# Input: Excel file with one row per protein and the columns
#   Group           protein accession (only needed to highlight proteins)
#   log2FoldChange  log2 fold change
#   Qvalue          adjusted p-value
# Decimal commas are accepted.
#
# Optional: a plain-text file with one accession per line; those proteins are
# drawn in dark grey on top of the rest.
#
# Usage: edit the settings below, then run  Rscript volcano_plot.R
# =============================================================================

for (pkg in c("readxl", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readxl)
library(ggplot2)

# --- Settings ----------------------------------------------------------------
input_file     <- "VolcanoPlot.xlsx"
highlight_file <- NULL            # e.g. "proteins_to_highlight.txt"
output_file    <- "volcano_plot.tiff"

fc_threshold <- 0.58              # |log2 fold change|
q_threshold  <- 0.05              # Q-value
x_limits     <- c(-5, 5)
y_limits     <- c(0, 10)

colours <- c("Up" = "#8e44ad", "Down" = "#00bfc4", "NS" = "grey80",
             "Highlight" = "#383d40")

# --- Data --------------------------------------------------------------------
d <- read_excel(input_file, col_types = "text")
d$Qvalue         <- as.numeric(gsub(",", ".", d$Qvalue))
d$log2FoldChange <- as.numeric(gsub(",", ".", d$log2FoldChange))
d$negLog10Q      <- -log10(d$Qvalue)

d$Significance <- with(d, ifelse(
  Qvalue < q_threshold & log2FoldChange >  fc_threshold, "Up",
  ifelse(Qvalue < q_threshold & log2FoldChange < -fc_threshold, "Down", "NS")
))

highlighted <- character(0)
if (!is.null(highlight_file)) highlighted <- trimws(readLines(highlight_file, warn = FALSE))
d$ColourGroup <- if ("Group" %in% names(d)) {
  ifelse(d$Group %in% highlighted, "Highlight", d$Significance)
} else d$Significance

message(sum(d$Significance == "Up"), " up, ", sum(d$Significance == "Down"),
        " down, ", sum(d$Significance == "NS"), " not significant")

# --- Plot --------------------------------------------------------------------
p <- ggplot() +
  geom_point(data = subset(d, ColourGroup != "Highlight"),
             aes(x = log2FoldChange, y = negLog10Q, colour = ColourGroup),
             alpha = 0.8, size = 1.8) +
  geom_point(data = subset(d, ColourGroup == "Highlight"),
             aes(x = log2FoldChange, y = negLog10Q),
             colour = colours[["Highlight"]], size = 2.2, alpha = 0.9) +
  scale_colour_manual(values = colours, breaks = c("Up", "Down", "NS")) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(q_threshold), linetype = "dashed") +
  labs(x = expression(Log[2]~Fold~Change),
       y = expression(-Log[10]~Q~value),
       colour = "Significance") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1.2),
        legend.position = "top")

tiff(output_file, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
print(p)
invisible(dev.off())
message("Saved ", output_file)
