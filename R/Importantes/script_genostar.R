library(GenoStaR)
library(readxl)
library(dplyr)
library(writexl)
ruta <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

# Leer archivos
df <- read_excel(file.path(ruta, "MATRIZ_FINAL_COMPLETA.xlsx"))

df <- df %>%
  mutate(across(everything(), ~ ifelse(is.na(.), " ", .)))


# SI PONGO ESPACIOS ME DA ERROR EN -CYP2D6_rs3892097
#df[is.na(df)] <- " " #Cambiar los NA por espacios para que funcione find_missing_data
#df[df == ""] <- " "
#df <- df %>% mutate(across(everything(), ~ ifelse(is.na(.), " ", .)))



df_sin_rs <- df %>%
  select(-CYP2D6_rs1065852)#, -CYP2D6_rs3892097)



lista_diplotipos <- assign_diplotype(df_sin_rs, c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6"), phased = FALSE, CYP1A2_name = "new")



diplotipos <- lista_diplotipos[[1]]

diplotipos <- diplotipos %>%
  dplyr::select(
    LabID.V2,
    CYP3A5_diplotype,
    CYP1A2_diplotype,
    CYP1A2_alternate_diplotype,
    CYP3A4_diplotype,
    CYP3A4_alternate_diplotype,
    CYP2C9_diplotype,
    CYP2B6_diplotype,
    CYP2B6_alternate_diplotype,
    CYP2C19_diplotype,
    CYP2D6_diplotype
  )

