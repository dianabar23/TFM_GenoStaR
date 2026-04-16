library(GenoStaR)
data(CYP2C9_data)




### PARA CYP3A4: sobre TODOS los individuos, luego se hace el genostar
#1. Cargar la tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/CYP3A5")
load("CYP3A5_Allele_def.rda")

#2. Cargar el csv con todos los individuos 
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/R")
load("matrix_geno_fixed_corregido.RData")
genotypes <- matrix_geno_fixed

df_CYP3A5 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A5_rs"))


#3. Filtrar para los SNPs en comun 
library(dplyr)
library(tidyr)

# 1. Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP3A5_Allele_Def)[-1]

# 2. Transformación y Filtrado Dinámico
df_resultado <- df_CYP3A5 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP3A5_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP3A5_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )




#4. Lanzar genostar solo paralos SNPs comunes (todos los individuos)
salida_CYP3A5_diplotipos <- assign_diplotype(df_resultado, 
                                             c("CYP3A5"),
                                             phased = FALSE)

df <- salida_CYP3A5_diplotipos[[1]] %>%
  select(LabID.V2, CYP3A5_diplotype)

df <- as.data.frame(df)

str(df$CYP3A5_diplotype)

salida_CYP3A5_metabolizador <- star_to_pheno(df, 
                                             c("CYP3A5"))


pie <- pie_chart(salida_CYP3A5_metabolizador, c("CYP3A5"))

