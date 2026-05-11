library(GenoStaR)
library(readxl)
library(dplyr)
library(writexl)
ruta <- "/Users/dianabarraso/Desktop/TFM_GenoStaR"

# Leer archivos
df <- read_excel(file.path(ruta, "MATRIZ_FINAL_COMPLETA.xlsx"))


lista_diplotipos <- assign_diplotype_diana(df, c(
  "CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6",
  "CYP2C19",
  "CYP2D6"), phased = FALSE
  , CYP1A2_name = "new"
  )

df_diplotipos <- lista_diplotipos[[1]]


### ACTIVITY CYP2C9 Y CYP2D6 
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



### PIE CHART 
genes <-c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6")

graficos <- pie_chart_diana(df_pheno, genes)


for (gen in genes) {
  
  # Construimos la ruta a la carpeta del gen (ej: .../TFM_GenoStaR/CYP2C9)
  carpeta_gen <- paste0(ruta, "/", gen)
  
  # Verificamos si la carpeta existe. Si no, la crea (por seguridad)
  if (!dir.exists(carpeta_gen)) {
    dir.create(carpeta_gen, recursive = TRUE)
  }
  
  # 3. Definimos el nombre y la ruta completa del archivo
  # Resultado: "/Users/.../TFM_GenoStaR/CYP2C9/CYP2C9_pie_chart.png"
  nombre_archivo <- paste0(carpeta_gen, "/", gen, "_pie_chart.png")
  
  # 4. Abrimos el dispositivo gráfico con la ruta completa
  png(filename = nombre_archivo, width = 800, height = 800, res = 120)
  
  # Ajustamos márgenes por si acaso para evitar el error anterior
  par(mar = c(4, 4, 4, 4))
  
  # 5. Llamamos a tu función
  pie_chart_diana(df_pheno, gen)
  
  # 6. CERRAMOS el archivo
  dev.off()
  
  message("Grafico guardado con éxito en: ", nombre_archivo)
}



### ALL GENO PHENO 
prueba_todos <- all_geno_pheno_diana(df, c(
  "CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6",
  "CYP2C19",
  "CYP2D6"), phased = FALSE
  , CYP1A2_name = "new"
)

df_final <- prueba_todos[[1]]
df_final <- as.data.frame(df_final)

# Guardar todo el dataframe de salida 
write_xlsx(
  df_final,
  path = file.path(ruta, "salida_genostar.xlsx")
)

# Guardar solo el dataframe con la salida 
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



### PIE CHART 
genes <-c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6")



for (gen in genes) {
  
  # Construimos la ruta a la carpeta del gen (ej: .../TFM_GenoStaR/CYP2C9)
  carpeta_gen <- paste0(ruta, "/", gen)
  
  # Verificamos si la carpeta existe. Si no, la crea (por seguridad)
  if (!dir.exists(carpeta_gen)) {
    dir.create(carpeta_gen, recursive = TRUE)
  }
  
  # 3. Definimos el nombre y la ruta completa del archivo
  # Resultado: "/Users/.../TFM_GenoStaR/CYP2C9/CYP2C9_pie_chart.png"
  nombre_archivo <- paste0(carpeta_gen, "/", gen, "_pie_chart.png")
  
  # 4. Abrimos el dispositivo gráfico con la ruta completa
  png(filename = nombre_archivo, width = 800, height = 800, res = 120)
  
  # Ajustamos márgenes por si acaso para evitar el error anterior
  par(mar = c(4, 4, 4, 4))
  
  # 5. Llamamos a tu función
  pie_chart_diana(df_final, gen)
  
  # 6. CERRAMOS el archivo
  dev.off()
  
  message("Grafico guardado con éxito en: ", nombre_archivo)
}

