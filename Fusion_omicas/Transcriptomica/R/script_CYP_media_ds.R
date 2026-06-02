library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

setwd("Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica")

gene_counts <- read_excel("CYP_gene_counts_normalizadas_nombres_transcriptomica_238_individuos.xlsx")


# =========================================================
# Calcular media y desviación típica por gen
# Data frame: gene_counts
# - Columna 1: CYP_name
# - Columna 2: EnsemblID
# - Columnas restantes: individuos
# =========================================================

# Seleccionar solo las columnas numéricas de individuos
expr_matrix <- gene_counts[, -(1:2)]

# Asegurarse de que todas las columnas son numéricas
expr_matrix <- as.data.frame(lapply(expr_matrix, as.numeric))

# Calcular media por gen (fila)
media_gen <- rowMeans(expr_matrix, na.rm = TRUE)

# Calcular desviación típica por gen
sd_gen <- apply(expr_matrix, 1, sd, na.rm = TRUE)

# Crear data frame final
resultados <- data.frame(
  CYP_name = gene_counts$CYP_name,
  EnsemblID = gene_counts$EnsemblID,
  Media = media_gen,
  Desviacion_tipica = sd_gen
)

# Ver resultados
print(resultados)


