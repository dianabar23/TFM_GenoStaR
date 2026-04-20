library(GenoStaR)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

#--------------------------------------------------
# 1. CONFIGURACIÓN
#--------------------------------------------------

base_dir <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

genes_config <- list(
  
  CYP1A2 = list(
    carpeta = "CYP1A2",
    file = "CYP1A2_Allele_def_2.rda",
    objeto = "CYP1A2_Allele_def_2"
  ),
  
  CYP3A4 = list(
    carpeta = "CYP3A4",
    file = "CYP3A4_Allele_def.rda",
    objeto = "CYP3A4_Allele_def"
  ),
  
  CYP3A5 = list(
    carpeta = "CYP3A5",
    file = "CYP3A5_Allele_def.rda",
    objeto = "CYP3A5_Allele_Def"
  ),
  
  CYP2C19 = list(
    carpeta = "CYP2C19",
    file = "CYP2C19_Allele_def.rda",
    objeto = "CYP2C19_Allele_def"
  ),
  
  CYP2C9 = list(
    carpeta = "CYP2C9",
    file = "CYP2C9_Allele_def.rda",
    objeto = "CYP2C9_Allele_def"
  ),
  
  CYP2B6 = list(
    carpeta = "CYP2B6",
    file = "CYP2B6_Allele_def.rda",
    objeto = "CYP2B6_Allele_def"
  ),
  
  CYP2D6 = list(
    carpeta = "CYP2D6",
    file = "CYP2D6_Allele_def.rda",
    objeto = "CYP2D6_Allele_def"
  )
)

#--------------------------------------------------
# 2. CARGAR GENOTIPOS
#--------------------------------------------------

ruta_genotipos <- file.path(base_dir, "matrix_geno_fixed_corregido.RData")

if (!file.exists(ruta_genotipos)) {
  stop(paste("❌ No se encuentra:", ruta_genotipos))
}

load(ruta_genotipos)
genotypes <- matrix_geno_fixed

#--------------------------------------------------
# 3. LISTA PARA MATRIZ FINAL GLOBAL
#--------------------------------------------------

lista_matrices <- list()

#--------------------------------------------------
# 4. FUNCIÓN DE PROCESAMIENTO
#--------------------------------------------------

procesar_gen <- function(gen, config) {
  
  cat("\n====================\nProcesando:", gen, "\n====================\n")
  
  # ruta completa del archivo .rda
  ruta_archivo <- file.path(base_dir, config$carpeta, config$file)
  
  if (!file.exists(ruta_archivo)) {
    stop(paste("❌ No se encuentra el archivo:", ruta_archivo))
  }
  
  # cargar referencia
  load(ruta_archivo)
  allele_def <- get(config$objeto)
  
  # seleccionar SNPs del gen
  df_gen <- genotypes %>% 
    select(LabID.V2, starts_with(paste0(gen, "_rs")))
  
  # SNPs de referencia
  cols_referencia <- colnames(allele_def)[-1]
  
  # pasar a formato largo y filtrar SNPs comunes
  df_long <- df_gen %>%
    pivot_longer(
      cols = starts_with(paste0(gen, "_")),
      names_to = "nombre_completo",
      values_to = "genotipo"
    ) %>%
    mutate(rs_limpio = sub(paste0(gen, "_"), "", nombre_completo)) %>%
    filter(rs_limpio %in% cols_referencia)
  
  # volver a formato ancho
  df_wide <- df_long %>%
    pivot_wider(
      id_cols = LabID.V2,
      names_from = nombre_completo,
      values_from = genotipo
    )
  
  return(df_wide)
}

#--------------------------------------------------
# 5. BUCLE PRINCIPAL
#--------------------------------------------------

for (gen in names(genes_config)) {
  
  df_gen <- procesar_gen(gen, genes_config[[gen]])
  
  lista_matrices[[gen]] <- df_gen
}

#--------------------------------------------------
# 6. MATRIZ FINAL GLOBAL
#--------------------------------------------------

df_final <- Reduce(function(x, y) {
  full_join(x, y, by = "LabID.V2")
}, lista_matrices)

#--------------------------------------------------
# 7. GUARDAR RESULTADO FINAL
#--------------------------------------------------

write_xlsx(
  df_final,
  file.path(base_dir, "MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx")
)

cat("\n✔ PROCESO COMPLETO TERMINADO\n")
cat("✔ Matriz final guardada como: MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx\n")
