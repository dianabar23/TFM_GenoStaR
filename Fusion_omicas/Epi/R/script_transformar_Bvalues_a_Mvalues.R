library(readxl)
library(dplyr)
library(writexl)

archivo <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGs_por_CYP_ventana100kb_con_betas.xlsx"

cyp_genes <- c("CYP2D6", "CYP2B6", "CYP1A2", "CYP3A4", "CYP3A5", "CYP2C9", "CYP2C19")

epi_list <- lapply(cyp_genes, function(gene) {
  read_excel(archivo, sheet = gene)
})

names(epi_list) <- cyp_genes




beta_to_m <- function(df) {
  
  df <- as.data.frame(df)
  
  # identificar columnas numéricas (muestras)
  idx <- sapply(df, is.numeric)
  
  # evitar dividir columnas de anotación
  df[idx] <- log2(df[idx] / (1 - df[idx]))
  
  return(df)
}

# aplicar a toda la lista
epi_list_M <- lapply(epi_list, beta_to_m)
names(epi_list_M) <- names(epi_list)



write_xlsx(
  epi_list_M,
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGs_por_CYP_ventana100kb_con_Mvalues.xlsx"
)
