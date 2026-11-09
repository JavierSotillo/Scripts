# =============================================================================
# Principal component analysis (PCA) of samples from a quantification table
#
# Input: tab-delimited file with proteins (or genes) in rows and samples in
# columns. The first column holds the identifiers; the remaining columns hold
# quantities for each sample. Decimal commas are accepted.
#
# Usage: edit the settings below, then run  Rscript pca_plot.R
# =============================================================================

for (pkg in c("readr", "ggplot2", "factoextra")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readr)
library(ggplot2)
library(factoextra)

# --- Settings ----------------------------------------------------------------
input_file  <- "PCA_Data.tsv"
output_file <- "PCA_plot.tiff"

# One label per sample column, in the same order as in the input file
group_labels <- factor(c("Control", "Control", "Control",
                         "Treated", "Treated", "Treated"))

# Normal-theory ellipses need at least 4 samples per group; with fewer,
# a convex hull is drawn around each group instead
ellipse_type <- if (min(table(group_labels)) >= 4) "norm" else "convex"

# --- Data --------------------------------------------------------------------
data <- read_delim(input_file, delim = "\t",
                   locale = locale(decimal_mark = ","), show_col_types = FALSE)
mat <- as.matrix(data[, -1])
rownames(mat) <- data[[1]]

if (ncol(mat) != length(group_labels)) {
  stop("The input has ", ncol(mat), " sample columns but ",
       length(group_labels), " group labels were given.")
}

# Scaled PCA needs complete rows with non-zero variance
keep <- complete.cases(mat) & apply(mat, 1, var, na.rm = TRUE) > 0
if (any(!keep)) message("Removed ", sum(!keep), " rows with missing values or zero variance")
mat <- mat[keep, ]

# --- PCA ---------------------------------------------------------------------
pca_result <- prcomp(t(mat), scale. = TRUE)
var_explained <- 100 * pca_result$sdev^2 / sum(pca_result$sdev^2)
message(sprintf("PC1: %.1f%%  PC2: %.1f%%", var_explained[1], var_explained[2]))

p <- fviz_pca_ind(pca_result,
                  geom.ind      = c("point", "text"),
                  col.ind       = group_labels,
                  addEllipses   = TRUE,
                  ellipse.type  = ellipse_type,
                  ellipse.level = 0.95,
                  palette       = "jco",
                  repel         = TRUE,
                  legend.title  = "Group",
                  title         = "")

tiff(output_file, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
print(p)
invisible(dev.off())
message("Saved ", output_file)
