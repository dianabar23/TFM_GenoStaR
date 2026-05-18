

#Cargar CYP1A2_resultados_final_diplotipos
#Cargar CYP3A4_resultados_final_diplotipos
#Cargar CYP3A5_resultados_final_diplotipos
#Cargar CYP2C19_resultados_final_diplotipos
#Cargar CYP2C9_resultados_final_diplotipos
#Cargar CYP2B6_resultados_final_diplotipos
#Cargar CYP2D6_resultados_final_diplotipos

library(dplyr)
library(ggplot2)
library(writexl)

# Ruta base
base_path <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

# Lista de genes
cyp_genes <- c("CYP1A2","CYP3A4","CYP3A5","CYP2C19","CYP2C9","CYP2B6","CYP2D6")

for (gene in cyp_genes) {
  
  obj_name <- paste0(gene, "_resultados_final_diplotipos")
  
  # Comprobar si existe el objeto
  if (!exists(obj_name)) {
    message(paste("No existe:", obj_name, "-> salto"))
    next
  }
  
  obj <- get(obj_name)
  
  # Crear ruta específica del CYP
  output_dir <- file.path(base_path, gene)
  
  # Crear carpeta si no existe
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extraer datos
  df <- obj[[1]] %>%
    select(LabID.V2, paste0(gene, "_diplotype"))
  
  colnames(df) <- c("LabID.V2", "Diplotype")
  
  # Guardar Excel
  write_xlsx(df, file.path(output_dir, paste0(gene, "_tabla_diplotipos.xlsx")))
  
  # Frecuencias
  freq <- df %>%
    count(Diplotype, sort = TRUE)
  
  # Gráfico
  p <- ggplot(freq, aes(x = reorder(Diplotype, n), y = n)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("Frecuencia de diplotipos", gene),
      x = "Diplotipo",
      y = "Frecuencia"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Guardar gráfico
  ggsave(
    filename = file.path(output_dir, paste0(gene, "_plot_frec_diplotipos.png")),
    plot = p,
    width = 8,
    height = 6
  )
  
  message(paste("Procesado:", gene, "-> guardado en", output_dir))
}






