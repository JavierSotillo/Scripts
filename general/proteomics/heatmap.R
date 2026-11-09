# =============================================================================
# Clustered heatmap with two sample annotations
#
# Input: Excel file without a header row, laid out as
#   row 1      sample names
#   row 2      first annotation for each sample (e.g. sample group)
#   row 3      second annotation for each sample (e.g. a clinical feature)
#   row 4 on   one protein per row
#   column 1   protein names (rows 4 onwards)
#
# Values are scaled by row (z-scores) and clustered with Manhattan distance
# and Ward's method.
#
# Usage: edit the settings below, then run  Rscript heatmap.R
# =============================================================================

for (pkg in c("readxl", "pheatmap")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readxl)
library(pheatmap)

# --- Settings ----------------------------------------------------------------
input_file  <- "Heatmap.xlsx"
output_file <- "Heatmap.tiff"

annotation_names <- c("Sample group", "Presence of granuloma")
annotation_colours <- list(
  "Sample group"          = c("Control" = "#440154", "SCC" = "#21918c", "UCC" = "#fde725"),
  "Presence of granuloma" = c("Yes" = "#FDBF6F", "No" = "#B3E5FC")
)

# --- Data --------------------------------------------------------------------
raw <- read_excel(input_file, col_names = FALSE, .name_repair = "minimal")

sample_names <- as.character(unlist(raw[1, -1]))
annotation_col <- data.frame(factor(as.character(unlist(raw[2, -1]))),
                             factor(as.character(unlist(raw[3, -1]))))
names(annotation_col) <- annotation_names
rownames(annotation_col) <- sample_names

values <- raw[-(1:3), -1]
data_matrix <- matrix(as.numeric(gsub(",", ".", as.matrix(values))),
                      nrow = nrow(values),
                      dimnames = list(raw[[1]][-(1:3)], sample_names))

# --- Plot --------------------------------------------------------------------
p <- pheatmap(data_matrix,
              scale = "row",
              annotation_col = annotation_col,
              annotation_colors = annotation_colours,
              clustering_distance_rows = "manhattan",
              clustering_distance_cols = "manhattan",
              clustering_method = "ward.D",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
              show_rownames = FALSE,
              show_colnames = FALSE,
              annotation_names_col_fontsize = 25,
              silent = TRUE)

tiff(output_file, width = 8, height = 5.6, units = "in", res = 300, compression = "lzw")
grid::grid.newpage()
grid::grid.draw(p$gtable)
invisible(dev.off())
message("Saved ", output_file)
