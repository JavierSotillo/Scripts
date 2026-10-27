# ============================================================
# Heatmap CON TODAS LAS PROTEÍNAS CUANTIFICADAS (sin filtro de significancia)
# Cultivo de Anisakis: 8 días x 3 réplicas
# ============================================================
install.packages("pheatmap")
install.packages("dplyr")
install.packages("readxl")
install.packages("dendextend")

library(pheatmap)
library(dplyr)
library(readxl)
library(dendextend)

# ------------------------------------------------------------
# 0. RUTA DEL ARCHIVO EXCEL
# ------------------------------------------------------------
excel_path <- "C:/Users/José Antonio/Documents/JULIA/COSAS 1ºPAPER/Scripts_R_HeatMap_PLSDA/Ani_ESPs_Report_Scripts_R.xlsx"

# ------------------------------------------------------------
# 1. CARGAR MATRIZ DE INTENSIDADES (formato ancho, log2)
# ------------------------------------------------------------
raw <- read_excel(excel_path, sheet = "Matrix", col_names = TRUE)
raw <- as.data.frame(raw, stringsAsFactors = FALSE)
colnames(raw)[1] <- "Samples"

#La fila "Group" contiene el nombre asociado a cada día de cultivo 
group_row     <- raw[raw$Samples == "Group", -1]
expr          <- raw[!raw$Samples != "Group", ]

rownames(expr) <- expr$Samples
expr$Samples   <- NULL
expr <- as.data.frame(lapply(expr, as.numeric), row.names = rownames(expr))

metadata <- data.frame(
  Sample    = colnames(expr),
  Group     = trimws(as.character(unlist(group_row))),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$Sample

print(table(metadata$Group))
#Debería salir cada día con 3 muestras

# ------------------------------------------------------------
# 2. SIN FILTRO DE SIGNIFICANCIA: usar todas las proteínas cuantificadas
# ------------------------------------------------------------
cat("Nº total de proteínas cuantificadas:", nrow(expr), "\n")

# Comprobación: la matriz ya viene imputada desde Spectronaut,
# así que no debería haber NA. Se verifica antes de continuar.
n_na <- sum(is.na(expr))
cat("Nº de valores NA en la matriz:", n_na, "\n")

if (n_na > 0) {
  stop("Hay NAs en la matriz pese a la imputación de Spectronaut — revisar antes de continuar.")
}

expr_complete <- expr

# ------------------------------------------------------------
# 3. CONSTRUIR MATRIZ Y ANOTACIONES (todas las proteínas, todos los días)
# ------------------------------------------------------------
mat_todos <- as.matrix(expr_complete)

annotation_todos <- metadata[colnames(mat_todos), "Group", drop = FALSE]
colnames(annotation_todos) <- "Culture day"

#Los días se ordenan tal y como aparecen en excel (de izq-drch)

annotation_colors_todos <- list(
  "Culture day" = c("Day 1" = "#FF9AA2", "Day 2" = "#FFB88C",
                    "Day 3" = "#FFE08A", "Day 4" = "#B5EAA2",
                    "Day 5" = "#8FE3E0", "Day 6" = "#8FBCFF",
                    "Day 7" = "#B7A6F5", "Day 8" = "#E8A6E8")
)

annotation_todos$`Culture day` <- factor(
  annotation_todos$`Culture day`,
  levels = names(annotation_colors_todos$`Culture day`)
)


cat("NAs en Culture day:", sum(is.na(annotation_todos$`Culture day`)), "\n")

# ------------------------------------------------------------
# 4. CLUSTERING DE COLUMNAS (Sin reordenación posterior pq no hay grupo Control )
# ------------------------------------------------------------
col_dist <- dist(scale(t(mat_todos)), method = "euclidean")
col_hc   <- hclust(col_dist, method = "ward.D2")

# ------------------------------------------------------------
# 5. HEATMAP (todas las proteínas, sin filtro de significancia)
# ------------------------------------------------------------
pheatmap(
  mat_todos,
  scale                    = "row",
  cluster_cols             = col_hc,
  clustering_method        = "ward.D2",
  annotation_col           = annotation_todos,
  annotation_colors        = annotation_colors_todos,
  annotation_names_col = FALSE,
  show_rownames            = FALSE,
  show_colnames            = FALSE,
  main                     = "",
  color = colorRampPalette(c("#5E3C99", "white", "#E66101"))(100),
  filename                 = file.path(dirname(excel_path), "heatmap_todas_proteinas.png"),
  width                    = 10,
  height                   = 8
)
