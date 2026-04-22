library(GenoStaR)
library(readxl)

### Cargar df filtrada 
df <- read_xlsx("MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx")
df[is.na(df)] <- " " # 2. Sustituyes los NA por un espacio en blanco " "

### Cargar df raw 
setwd("/Users/dianabarraso/Desktop/TFM_GenoStaR")
load("matrix_geno_fixed_espacios.RData")
df_raw <- matrix_geno_fixed



### FUNCION: find_missing_data()
missing_data <- find_missing_data(df)
missing_data_raw <- find_missing_data(df_raw)


### FUNCION: fill_empty_cells()
setwd("/Users/dianabarraso/Desktop/TFM_GenoStaR/CYP2C19")
load("CYP2C19_Allele_def.rda")
CYP2C19_allele_def <- CYP2C19_Allele_def
CYP2C19_allele_def <- CYP2C19_allele_def[-1, ]
CYP2C19_allele_def_COMPLETA <- fill_empty_cells(CYP2C19_allele_def)
## Se rellena 



### FUNCION: find_diplotype(data, genes, CYP1A2_name = )
diplotipos_raw <- find_diplotype(df_raw, c("CYP1A2","CYP3A4", "CYP3A5", "CYP2B6", "CYP2C19","CYP2C9"), CYP1A2_name = "new")
diplotipos <- find_diplotype(df, c("CYP1A2","CYP3A4", "CYP3A5", "CYP2B6", "CYP2C19","CYP2C9"), CYP1A2_name = "new")


diplotipos_df_raw <- diplotipos_raw[[1]]
diplotipos_df <- diplotipos[[1]]

no_match <- filter_no_match(diplotipos_df, "CYP1A2")
