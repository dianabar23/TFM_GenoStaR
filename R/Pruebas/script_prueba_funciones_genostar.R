library(GenoStaR)
library(readxl)
library(dplyr)
library(purrr)


setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")


df_final <- read_xlsx("MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx")



#SACAR DIPLOTIPOS 
salida_diplotipos <- assign_diplotype(df_final, 
                                             c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19"),
                                             phased = FALSE, CYP1A2_name = "new")

genes <- c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19")

df_diplotipos_final <- map(genes, function(gen) {
  
  salida_diplotipos[[1]] %>%
    select(LabID.V2, any_of(paste0(gen, "_diplotype")))
  
}) %>%
  reduce(full_join, by = "LabID.V2")

df_diplotipos_final <- as.data.frame(df_diplotipos_final)

df_diplotipos_final <- df_diplotipos_final %>%
  mutate(across(-LabID.V2, as.character))


# SACAR ACTIVITY
genes_activity <- c("CYP2C9") #METER CYP2D6
salida_activity <- map(genes_activity, function(gen) {
  
  cat("Procesando:", gen, "\n")
  
  df <- star_to_activity(df_diplotipos_final, gen)
  
  # Nos quedamos solo con la columna relevante
  df %>%
    select(LabID.V2, contains(gen))
})

# Unir todos los resultados
df_activity_final <- reduce(salida_activity, full_join, by = "LabID.V2")




# SACAR METABOLIZADOR 
genes_metabolizador <- c("CYP3A5","CYP2C9","CYP2B6","CYP2C19") #METER CYP2D6
salida_metabolizador <- map(genes_metabolizador, function(gen) {
  
  cat("Procesando:", gen, "\n")
  
  df <- star_to_pheno(df_diplotipos_final, gen)
  
  # Nos quedamos solo con la columna relevante
  df %>%
    select(LabID.V2, contains(gen))
})

# Unir todos los resultados
df_metabolizador_final <- reduce(salida_metabolizador, full_join, by = "LabID.V2")




# PIE CHART

salida_metabolizador <- star_to_pheno(df_diplotipos_final, 
                                      c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19"))


pie <- pie_chart(salida_CYP3A5_metabolizador, c("CYP3A5"))


all_geno_pheno(df_resultado, 
                c("CYP3A5"),
                phased = FALSE)
