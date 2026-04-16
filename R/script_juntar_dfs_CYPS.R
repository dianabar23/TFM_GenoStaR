library(GenoStaR)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

#--------------------------------------------------
# 1. CONFIGURACIÓN
#--------------------------------------------------

setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")

genes_config <- list(
  CYP1A2 = list(file = "CYP1A2_Allele_def_2.rda", objeto = "CYP1A2_Allele_def_2"),
  CYP3A4 = list(file = "CYP3A4_Allele_def.rda", objeto = "CYP3A4_Allele_def"),
  CYP3A5 = list(file = "CYP3A5_Allele_def.rda", objeto = "CYP3A5_Allele_Def"),
  CYP2C19 = list(file = "CYP2C19_Allele_def.rda", objeto = "CYP2C19_Allele_def"),
  CYP2C9 = list(file = "CYP2C9_Allele_def.rda", objeto = "CYP2C9_Allele_def"),
  CYP2B6 = list(file = "CYP2B6_Allele_def.rda", objeto = "CYP2B6_Allele_def"),
  CYP2D6 = list(file = "CYP2D6_Allele_def.rda", objeto = "CYP2D6_Allele_def")
)

#--------------------------------------------------
# 2. CARGAR GENOTIPOS
#--------------------------------------------------

load("matrix_geno_fixed_corregido.RData")
genotypes <- matrix_geno_fixed

#--------------------------------------------------
# 3. LISTA PARA MATRIZ FINAL GLOBAL
#--------------------------------------------------

lista_matrices <- list()

#--------------------------------------------------
# 4. FUNCIÓN DE PREPARACIÓN (SIN GENOSTAR)
#--------------------------------------------------

procesar_gen <- function(gen, config) {
  
  cat("\n====================\nProcesando:", gen, "\n====================\n")
  
  # cargar referencia
  load(config$file)
  allele_def <- get(config$objeto)
  
  # subset genotipos
  df_gen <- genotypes %>% 
    select(LabID.V2, starts_with(paste0(gen, "_rs")))
  
  # SNPs de referencia
  cols_referencia <- colnames(allele_def)[-1]
  
  #---------------------------------------
  # FILTRADO A SNPs COMUNES
  #---------------------------------------
  
  df_long <- df_gen %>%
    pivot_longer(
      cols = starts_with(paste0(gen, "_")),
      names_to = "nombre_completo",
      values_to = "genotipo"
    ) %>%
    mutate(rs_limpio = sub(paste0(gen, "_"), "", nombre_completo)) %>%
    filter(rs_limpio %in% cols_referencia)
  
  df_wide <- df_long %>%
    pivot_wider(
      id_cols = LabID.V2,
      names_from = nombre_completo,
      values_from = genotipo
    )
  
  #---------------------------------------
  # GUARDAR EXCEL POR GEN
  #---------------------------------------
  
  write_xlsx(
    df_wide,
    paste0(gen, "_genotipos_filtrados_comunes.xlsx")
  )
  
  # devolver para matriz global
  return(df_wide)
}

#--------------------------------------------------
# 5. BUCLE PRINCIPAL
#--------------------------------------------------

for (gen in names(genes_config)) {
  
  df_gen <- procesar_gen(gen, genes_config[[gen]])
  
  # guardar en lista
  lista_matrices[[gen]] <- df_gen
}

#--------------------------------------------------
# 6. MATRIZ FINAL GLOBAL (TODOS LOS GENES)
#--------------------------------------------------

# merge progresivo por LabID.V2
df_final <- Reduce(function(x, y) {
  full_join(x, y, by = "LabID.V2")
}, lista_matrices)

#--------------------------------------------------
# 7. GUARDAR MATRIZ GLOBAL
#--------------------------------------------------

write_xlsx(
  df_final,
  "MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx"
)

cat("\n✔ PROCESO COMPLETO TERMINADO\n")
cat("✔ Matriz final guardada como: MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx\n")