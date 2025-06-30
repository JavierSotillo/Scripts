# Load necessary libraries
library(ggplot2)
library(readr)
library(factoextra)

#Set working directory
setwd("C:/Users/javier.sotillo/Desktop")

# Set file path
file_path <- "Heatmap_Data.tsv"

# Step 1: Load the data
data <- read_delim("Heatmap_Data.tsv", 
                   delim = "\t", 
                   locale = locale(decimal_mark = ","))

# Step 2: Set protein names as row names, remove first column
rownames(data) <- data[[1]]
data <- data[ , -1]

# Step 3: Transpose: now rows = samples, columns = proteins
data_t <- t(data)

# Step 4: Manually define group labels
# Example: 5 replicates each for A, B, C
group_labels <- factor(c("WT","IL10KO","WT","WT","IL10KO","IL10KO","IL10KO","IL10KO","WT","WT","IL10KO_AHCC","IL10KO_AHCC","IL10KO_AHCC","IL10KO_AHCC","IL10KO_AHCC"))

# Step 5: PCA
pca_result <- prcomp(data_t, scale. = TRUE)

# Step 6: Plot with ellipses
p <- fviz_pca_ind(pca_result,
                  geom.ind = "point",
                  col.ind = group_labels,
                  addEllipses = TRUE,
                  ellipse.level = 0.95,
                  palette = "jco",
                  repel = TRUE,
                  legend.title = "Group")

# Step 7: Save and show plot
ggsave("PCA_manual_groups.png", plot = p, width = 8, height = 6, dpi = 300)
print(p)