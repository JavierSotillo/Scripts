### PCA PLOT desde Excel: 8 días de cultivo x 3 réplicas ####

library(readxl)
library(ggplot2)

# 1. Cargar archivo Excel
#    >>> CAMBIA esta ruta por la de tu Excel (usa "/" en lugar de "\") <<<
file_path <- "C:/Users/José Antonio/Documents/JULIA/COSAS 1ºPAPER/Scripts_R_HeatMap_PLSDA/Ani_ESPs_Report_Scripts_R.xlsx"
data_raw <- read_excel(file_path, col_names = FALSE)  # leer sin nombres de columna automáticos

# 2. Localizar la fila "Group" y extraer las etiquetas de grupo (sin la primera columna)
fila_group <- which(tolower(trimws(as.character(data_raw[[1]]))) == "group")
if (length(fila_group) != 1) {
  stop("❗ No encuentro una (y solo una) fila llamada 'Group' en la primera columna.")
}

group_labels <- as.character(unlist(data_raw[fila_group, -1]))
group_labels <- factor(group_labels, levels = unique(group_labels))  # mantiene el orden Day 1 ... Day 8

# Nombres de las muestras (primera fila = cabecera)
sample_names <- as.character(unlist(data_raw[1, -1]))

# 3. Eliminar la cabecera (fila 1) y la fila de grupos
data <- data_raw[-c(1, fila_group), ]

# 4. Convertir a data.frame y establecer nombres de fila (proteínas)
data <- as.data.frame(data)
rownames(data) <- data[[1]]  # primera columna = nombres de proteína
data <- data[, -1]

# 5. Convertir comas a puntos y a numérico
data[] <- lapply(data, function(x) as.numeric(gsub(",", ".", as.character(x))))

# 6. Transponer: filas = muestras, columnas = proteínas
data_t <- t(data)

# 7. Eliminar columnas (proteínas) con NA o varianza cero
data_t <- data_t[, colSums(is.na(data_t)) == 0, drop = FALSE]
data_t <- data_t[, apply(data_t, 2, var, na.rm = TRUE) != 0, drop = FALSE]

# 8. Asegurarse que group_labels tiene el mismo orden que las filas transpuestas
if (nrow(data_t) != length(group_labels)) {
  stop("❗ Número de etiquetas de grupo no coincide con número de muestras tras transposición.")
}

# 9. Ejecutar PCA
if (ncol(data_t) > 0) {
  pca_result <- prcomp(data_t, scale. = TRUE)

  # Porcentaje de varianza explicada
  var_exp <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

  # Datos del PCA
  pca_coords <- as.data.frame(pca_result$x[, 1:2])
  colnames(pca_coords)[1:2] <- c("PCA1", "PCA2")
  pca_coords$Group  <- group_labels
  pca_coords$Sample <- sample_names

  # Elipses calculadas a mano para cada día (stat_ellipse no funciona con solo 3 réplicas).
  # Cada elipse está centrada en la media de las 3 réplicas y su forma sale de su covarianza.
  nivel_elipse <- 0.90   # cambia a 0.90, 0.68... si quieres una elipse más pequeña

  calcular_elipse <- function(d, nivel = nivel_elipse, n = 100) {
    centro <- colMeans(d[, c("PCA1", "PCA2")])
    covarianza <- cov(d[, c("PCA1", "PCA2")])
    radio <- sqrt(qchisq(nivel, df = 2))
    ev <- eigen(covarianza)
    A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))
    angulo <- seq(0, 2 * pi, length.out = n)
    puntos <- cbind(cos(angulo), sin(angulo)) %*% t(A)
    puntos <- sweep(radio * puntos, 2, centro, "+")
    data.frame(PCA1 = puntos[, 1], PCA2 = puntos[, 2], Group = d$Group[1])
  }

  elipses <- do.call(rbind, lapply(split(pca_coords, pca_coords$Group), calcular_elipse))

  # 8 colores tipo arcoíris en tono pastel (uno por día; puedes cambiarlos)
  colores_dias <- c("#F28B82",  # Day 1 - rojo pastel
                    "#F8B26A",  # Day 2 - naranja pastel
                    "#F3D55B",  # Day 3 - amarillo pastel
                    "#9BD17F",  # Day 4 - verde pastel
                    "#6CCDBE",  # Day 5 - turquesa pastel
                    "#7FB6EA",  # Day 6 - azul pastel
                    "#9E94E0",  # Day 7 - índigo pastel
                    "#D58FD0")  # Day 8 - violeta pastel
  names(colores_dias) <- levels(group_labels)

  # Plot con ggplot2
  p <- ggplot(pca_coords, aes(x = PCA1, y = PCA2, color = Group)) +
    # Elipses por día (calculadas arriba; stat_ellipse necesita al menos 4 puntos por grupo)
    geom_polygon(data = elipses, aes(fill = Group), alpha = 0.2, linewidth = 0.8) +
    geom_point(shape = 16, size = 6) +                              # círculos iguales
    scale_color_manual(values = colores_dias) +
    scale_fill_manual(values = colores_dias) +
    labs(
      x = paste0("PCA1 (", var_exp[1], "%)"),
      y = paste0("PCA2 (", var_exp[2], "%)"),
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", size = 16),   # títulos de ejes en negrita y grandes
      axis.text = element_text(size = 14),                   # números de los ejes más grandes
      legend.title = element_text(face = "bold", size = 14), # título de la leyenda en negrita
      legend.text = element_text(size = 13)                  # texto de la leyenda más grande
    )

  # Guardar la figura en la misma carpeta que el Excel
  ggsave(file.path(dirname(file_path), "PCA_8dias.png"),
         plot = p, width = 10, height = 7, dpi = 300)
  print(p)
}
