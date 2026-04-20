library(GenoStaR)
library(readxl)
library(dplyr)
library(purrr)

setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")

#---------------------------------------
# CARGAR DATOS
#---------------------------------------

df_final <- read_xlsx("MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx")

#---------------------------------------
# DIPLOTIPOS
#---------------------------------------

salida_diplotipos <- assign_diplotype(
  df_final, 
  c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19"),
  phased = FALSE,
  CYP1A2_name = "new"
)

genes <- c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19")

df_diplotipos <- map(genes, function(gen) {
  
  salida_diplotipos[[1]] %>%
    select(LabID.V2, any_of(paste0(gen, "_diplotype")))
  
}) %>%
  reduce(full_join, by = "LabID.V2") %>%
  mutate(across(-LabID.V2, as.character))

#---------------------------------------
# FUNCIÓN
#---------------------------------------

procesar_gen_completo <- function(gen, df_diplotipos) {
  
  cat("\nProcesando:", gen, "\n")
  
  col_diplo <- paste0(gen, "_diplotype")
  
  if (!(col_diplo %in% colnames(df_diplotipos))) {
    cat("⏭️ No hay diplotipo para", gen, "\n")
    return(NULL)
  }
  
  # DIPLOTIPO
  df_base <- df_diplotipos %>%
    select(LabID.V2, all_of(col_diplo))
  
  # ACTIVITY
  df_activity <- tryCatch({
    star_to_activity(df_diplotipos, gen) %>%
      select(LabID.V2, contains(gen))
  }, error = function(e) NULL)
  
  # PHENOTYPE + PIE CHART
  df_pheno <- tryCatch({
    
    df_tmp <- star_to_pheno(df_diplotipos, gen)
    
    df_tmp_sel <- df_tmp %>%
      select(LabID.V2, contains(gen))
    
    col_pheno <- grep(gen, colnames(df_tmp_sel), value = TRUE)
    
    if (length(col_pheno) > 0 &&
        any(!is.na(df_tmp_sel[[col_pheno[1]]]))) {
      
      pie_chart(df_tmp_sel, c(gen))
    }
    
    df_tmp_sel
    
  }, error = function(e) NULL)
  
  # UNIR
  lista <- list(df_base, df_activity, df_pheno) %>%
    compact()
  
  reduce(lista, full_join, by = "LabID.V2")
}

#---------------------------------------
# PIPELINE FINAL
#---------------------------------------

df_final_completo <- map(
  genes,
  ~procesar_gen_completo(.x, df_diplotipos)
) %>%
  compact() %>%
  reduce(full_join, by = "LabID.V2")
