# --- Load packages ---
install.packages(c("readxl", "pheatmap"))
library(readxl)
library(pheatmap)

# --- Read the Excel file ---
# Assumes: 
#   - Row 1 = sample names
#   - Row 2 = first group labels
#   - Row 3 = second group labels
#   - Column 1 = protein names
raw <- read_excel("C:/Users/javier.sotillo/OneDrive - Instituto de Salud Carlos III/Papers/Paper_FFPE_Schisto/Figures/Heatmap/FFPE_Heatmap.xlsx", col_names = FALSE)

# --- Extract metadata and data ---
sample_names <- as.character(unlist(raw[1, -1]))   # first row, excluding first column
group1 <- as.character(unlist(raw[2, -1]))         # second row, first group
group2 <- as.character(unlist(raw[3, -1]))         # third row, second group
data_matrix <- as.data.frame(raw[-c(1,2,3), -1])   # data from row 4 onward

# --- Assign proper row/column names ---
rownames(data_matrix) <- raw[-c(1,2,3), 1, drop = TRUE]  # protein names
colnames(data_matrix) <- sample_names

# --- Convert to numeric matrix ---
data_matrix <- as.matrix(sapply(data_matrix, as.numeric))
rownames(data_matrix) <- raw[-c(1,2,3), 1, drop = TRUE]

# --- Create column annotation for groups ---
annotation_col <- data.frame(
  "Sample group" = factor(group1),
  "Presence of granuloma" = factor(group2),
  check.names = FALSE  # <- esto conserva los espacios
)
rownames(annotation_col) <- sample_names

# --- Define colors for groups ---
group_colors <- list(
  "Sample group" = c("Control" = "#440154", "SCC" = "#21918c", "UCC" = "#fde725"),
  "Presence of granuloma" = c("Yes" = "#FDBF6F", "No" = "#B3E5FC")
)

# --- Plot the heatmap ---
p <- pheatmap(data_matrix,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = group_colors,
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         clustering_method = "ward.D",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE,
         silent = TRUE,
         annotation_names_col_fontsize = 25  # tamaño de fuente aumentado sin usar gpar
)

print(p)


tiff("C:/Users/javier.sotillo/OneDrive - Instituto de Salud Carlos III/Papers/Paper_FFPE_Schisto/Figures/Heatmap/Heatmap.tiff", width = 8, height = 5.6, units = "in", res = 300, compression = "lzw")
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()
