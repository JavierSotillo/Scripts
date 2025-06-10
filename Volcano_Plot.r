# Instala y carga paquetes necesarios
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(readxl)
library(ggplot2)

# Definir la ruta del archivo
ruta_archivo <- "C:/Users/javier.sotillo/Desktop/VolcanoPlot.xlsx"  


# Asumimos que el archivo tiene columnas llamadas 'log2FoldChange' y 'Qvalue'
datos <- read_excel(ruta_archivo)

# Crear columna -log10(Qvalue)
datos$negLog10Qvalue <- -log10(datos$Qvalue)

# Umbrales de significancia
fc_threshold <- 0.58
Q_threshold <- 0.05

# Clasificación para colorear
datos$Significance <- with(datos, ifelse(
  Qvalue < Q_threshold & log2FoldChange > fc_threshold, "Up",
  ifelse(Qvalue < Q_threshold & log2FoldChange < -fc_threshold, "Down", "NS")
))

# Mapa de colores: morado (up), cian (down), gris (no significativo)
colores <- c("Up" = "#8e44ad", "Down" = "#00bfc4", "NS" = "grey80")

# Crear Volcano Plot
volcano_plot <- ggplot(datos, aes(x = log2FoldChange, y = negLog10Qvalue)) +
  geom_point(aes(color = Significance), alpha = 0.8, size = 1.8) +
  scale_color_manual(values = colores) +
  scale_x_continuous(limits = c(-5, 5), expand = c(0, 0)) +
scale_y_continuous(limits = c(0.5, 10), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Q~value),
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(Q_threshold), linetype = "dashed", color = "black") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 1.2),  # líneas de ejes negras y más gruesas
    legend.position = "top"
  )

# Exportar Volcano plot a TIFF
tiff("C:/Users/javier.sotillo/Desktop/volcano_plot.tiff", width = 8, height = 6, units = "in", res = 300)
print(volcano_plot)
dev.off()
