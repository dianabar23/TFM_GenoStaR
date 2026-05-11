library(GenoStaR)
library(readxl)
library(dplyr)
library(writexl)


#ruta <- "/Users/dianabarraso/Desktop/TFM_GenoStaR"
ruta <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR"

# Leer archivos
df <- read_excel(file.path(ruta, "MATRIZ_FINAL_COMPLETA.xlsx"))
df <- as.data.frame(df)


### FUNCION: ALL GENO PHENO 
df_todo <- all_geno_pheno_diana(df, c(
  "CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6",
  "CYP2C19",
  "CYP2D6"), phased = FALSE
  , CYP1A2_name = "new"
)


df_final <- df_todo[[1]]
df_final <- as.data.frame(df_final)

# Guardar todo el dataframe de salida 
write_xlsx(
  df_final,
  path = file.path(ruta, "salida_genostar.xlsx")
)

# Guardar solo el dataframe de salida con las columnas que nos interesan (filtrado)
df_final_filtrado <- df_final[, c(
  "LabID.V2",
  "CYP3A5_diplotype",
  "CYP3A5_AS_1",
  "CYP3A5_Metabolizer_Status",
  "CYP3A5_Comment",
  
  "CYP1A2_diplotype",
  "CYP1A2_alternate_diplotype",
  "CYP1A2_AS_1",
  "CYP1A2_Metabolizer_Status",
  "CYP1A2_Comment",
  
  "CYP3A4_diplotype",
  "CYP3A4_alternate_diplotype",
  "CYP3A4_AS_1",
  "CYP3A4_Metabolizer_Status",
  "CYP3A4_Comment",
  
  "CYP2C9_diplotype",
  "CYP2C9_AS_1",
  "CYP2C9_Metabolizer_Status",
  "CYP2C9_Comment",
  
  "CYP2B6_diplotype",
  "CYP2B6_alternate_diplotype",
  "CYP2B6_AS_1",
  "CYP2B6_Metabolizer_Status",
  "CYP2B6_Comment",
  
  "CYP2C19_diplotype",
  "CYP2C19_AS_1",
  "CYP2C19_Metabolizer_Status",
  "CYP2C19_Comment",
  
  "CYP2D6_diplotype",
  "CYP2D6_alternate_diplotype",
  "CYP2D6_AS_1",
  "CYP2D6_AS_2",
  "CYP2D6_Metabolizer_Status",
  "CYP2D6_Metabolizer_Status_2",
  "CYP2D6_Comment"
)]


write_xlsx(
  df_final_filtrado,
  path = file.path(ruta, "salida_genostar_filtrado.xlsx")
)



### FUNCION: PIE CHART
# Se necesita bucle porque va de gen en gen 
genes <-c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6")

for (gen in genes) {
  
  phenotype_col <- paste0(gen, "_Metabolizer_Status")
  
  # 1. Comprobar que la columna existe
  if (!phenotype_col %in% names(df_final)) {
    message("Columna no existe para: ", gen)
    next
  }
  
  # 2. Obtener datos
  datos <- df_final[[phenotype_col]]
  
  # 3. Quitar NA
  datos <- datos[!is.na(datos)]
  
  # 4. Comprobar si hay datos válidos
  if (length(datos) == 0) {
    message("Sin datos para: ", gen)
    next
  }
  
  # 5. Comprobar que no todo es 0 (por seguridad)
  if (sum(table(datos)) == 0) {
    message("Datos inválidos para: ", gen)
    next
  }
  
  # --- A partir de aquí ya es seguro ---
  
  carpeta_gen <- paste0(ruta, "/", gen)
  
  if (!dir.exists(carpeta_gen)) {
    dir.create(carpeta_gen, recursive = TRUE)
  }
  
  nombre_archivo <- paste0(carpeta_gen, "/", gen, "_pie_chart.png")
  
  png(filename = nombre_archivo, width = 1800, height = 1600, res = 200)
  
  par(mar = c(8, 8, 6, 8), cex = 0.7)
  
  pie_chart_diana(df_final, gen)
  
  dev.off()
  
  message("Grafico guardado: ", nombre_archivo)
}




#### AUNQUE LA FUNCION GENERAL FUNCIONA SE VAN PROBANDO LAS OTRAS FUNCIONES UNA A UNA 

### FUNCION: ASSIGN_DIPLOTYPE 

lista_diplotipos <- assign_diplotype_diana(df, c(
  "CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6",
  "CYP2C19",
  "CYP2D6"), phased = FALSE
  , CYP1A2_name = "new"
  )

df_diplotipos <- lista_diplotipos[[1]]



### FUNCION: STAR TO ACTIVITY para CYP2C9 Y CYP2D6 
# 1. Forzamos a que sea un data.frame puro (esto quita el comportamiento de Tibble)
df_diplotipos <- as.data.frame(df_diplotipos)

df_activity <- star_to_activity_diana(df_diplotipos, "CYP2C9")
df_activity <- star_to_activity_diana(df_activity, "CYP2D6")


### PHENO CYP2C19, CYP3A5, CYP2B6,  CYP2C9 Y CYP2D6 
# 1. Forzamos a que sea un data.frame puro (esto quita el comportamiento de Tibble)
df_activity <- as.data.frame(df_activity)
df_pheno <- star_to_pheno_diana(df_activity, "CYP2D6")
df_pheno <- star_to_pheno_diana(df_pheno, "CYP2B6")
df_pheno <- star_to_pheno_diana(df_pheno, "CYP2C19")
df_pheno <- star_to_pheno_diana(df_pheno, "CYP2C9")
df_pheno <- star_to_pheno_diana(df_pheno, "CYP3A5")


