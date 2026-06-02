### LOAD LIBRARIES
suppressMessages({
  library(readxl)
  library(ggplot2)
  library(writexl)
  library(car)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

##################
# Modelo gene_counts ~ metabolizador para cada CYP 
##################

### 1. CARGAR TIPO METABOLIZADOR DE LOS 238 INDIVIDUOS 
metabolizador <- read_excel("Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Farmacogenomica_genostar/CYP_metabolizador_nombres_genostar_660_individuos.xlsx")

metabolizador <- metabolizador %>%
  column_to_rownames("LabID.V2")

### 2. CARGAR GENE COUNTS
gene_counts <- read_excel(
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_transcriptomica_660_individuos.xlsx"
)

gene_counts <- gene_counts[, -2]

gene_counts <- gene_counts %>%
  column_to_rownames("CYP_name")

### 3. REORDENAR AMBAS VARIABLES PARA MISMO ORDEN DE MUESTRAS
common_samples <- intersect(colnames(gene_counts), rownames(metabolizador))

gene_counts <- gene_counts[, common_samples]
metabolizador <- metabolizador[common_samples, , drop = FALSE]

stopifnot(all(colnames(gene_counts) == rownames(metabolizador)))

### 4. PREPARAR GENE COUNTS 
gene_counts <- gene_counts[!rownames(gene_counts) %in% c("CYP1A2", "CYP3A4"), ]
gene_counts_scaled <- as.data.frame(t(scale(t(gene_counts))))

### 5. PREPARAR METABOLIZADORES (quitar NO e INDETERMINATE)
invalid <- c("NO", "Indeterminate")

metabolizador_clean <- metabolizador %>%
  mutate(across(everything(), ~ ifelse(.x %in% invalid, NA, .x)))

### 6. MODELOS LINEALES (UNO POR CYP)

cyps <- c("CYP2D6", "CYP2B6", "CYP3A5", "CYP2C9", "CYP2C19")

resultados <- list()

for (cyp in cyps) {
  
  df <- data.frame(
    expr = as.numeric(gene_counts_scaled[cyp, ]),
    metabolizador = metabolizador_clean[[cyp]]
  )
  
  # quitar NA
  df <- df[!is.na(df$metabolizador) & !is.na(df$expr), ]
  
  # ordenar factor biológico
  df$metabolizador <- factor(
    df$metabolizador,
    levels = c("Poor", "Intermediate", "Normal", "Rapid", "Ultrarapid")
  )
  
  # modelo lineal clásico
  fit <- lm(expr ~ metabolizador, data = df)
  
  # ANOVA + resumen
  anova_res <- anova(fit)
  coef_res <- summary(fit)
  
  resultados[[cyp]] <- list(
    model = fit,
    anova = anova_res,
    summary = coef_res
  )
}

### 7. EJEMPLO RESULTADO CYP2D6
resultados[["CYP2D6"]]
resultados[["CYP3A5"]]
resultados[["CYP2B6"]]
resultados[["CYP2C9"]]
resultados[["CYP2C19"]]

