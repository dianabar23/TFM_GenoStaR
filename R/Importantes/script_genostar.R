library(GenoStaR)
library(readxl)
library(dplyr)
library(writexl)
ruta <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

# Leer archivos
df <- read_excel(file.path(ruta, "MATRIZ_FINAL_COMPLETA.xlsx"))

#df <- df %>% mutate(across(everything(), ~ ifelse(is.na(.), " ", .)))



lista_diplotipos <- assign_diplotype_diana(df, c(
  "CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6",
  "CYP2C19",
  "CYP2D6"), phased = FALSE
  , CYP1A2_name = "new"
  )

### LLEGUE HASTA AQUI 

df_diplotipos <- lista_diplotipos[[1]]

df_diplotipos <- df_diplotipos %>%
  mutate(across(-LabID.V2, as.character))

diplotipos <- df_diplotipos %>%
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

diplotipos <- diplotipos %>%
  mutate(across(-LabID.V2, ~ na_if(., "NA")))

diplotipos <- mutate(diplotipos, across(-LabID.V2, as.character))


write_xlsx(diplotipos, "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/diplotipos.xlsx")


diplotipos <- read_xlsx("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/diplotipos.xlsx")
diplotipos[] <- lapply(diplotipos, function(x) {
  as.character(trimws(x))
})

genes <- c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6")

CYP2C9_activity <- star_to_pheno(lista_diplotipos, "CYP2C19")

                             