
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
library(scales)

# Ruta base
base_path <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

# Genes
cyp_genes <- c("CYP1A2","CYP3A4","CYP3A5","CYP2C19","CYP2C9","CYP2B6","CYP2D6")

for (gene in cyp_genes) {
  
  obj_name <- paste0(gene, "_resultados_final_diplotipos")
  
  # Si no existe, saltar
  if (!exists(obj_name)) {
    message(paste("No existe:", obj_name, "-> salto"))
    next
  }
  
  obj <- get(obj_name)
  
  # Carpeta salida
  output_dir <- file.path(base_path, gene)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Datos
  df <- obj[[1]] %>%
    select(LabID.V2, paste0(gene, "_diplotype"))
  
  colnames(df) <- c("LabID.V2", "Diplotype")
  
  # Guardar Excel
  write_xlsx(df, file.path(output_dir, paste0(gene, "_tabla_diplotipos.xlsx")))
  
  # Frecuencias
  # Frecuencias (NA como categoría normal)
  freq <- df %>%
    count(Diplotype) %>%
    mutate(
      prop = n / sum(n),
      categoria = case_when(
        Diplotype == "NA" ~ "NA",
        Diplotype == "*1/*1" ~ "Normal (*1/*1)",
        TRUE ~ "Otros"
      )
    ) %>%
    arrange(desc(prop))
  
  # Gráfico
  p <- ggplot(freq, aes(x = reorder(Diplotype, prop), y = prop, fill = categoria)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    geom_text(aes(label = ifelse(is.na(Diplotype), "NA", percent(prop, accuracy = 0.1))),
              hjust = -0.1, size = 3) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(
      values = c(
        "NA" = "#D73027",               # rojo
        "Normal (*1/*1)" = "#1A9850",   # verde
        "Otros" = "#4575B4"             # azul
      )
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Distribución de diplotipos", gene),
      subtitle = "Proporciones de diplotipos con NA y *1/*1 destacados",
      x = "Diplotipo",
      y = "Proporción",
      fill = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    expand_limits(y = max(freq$prop, na.rm = TRUE) * 1.15)
  
  # Guardar gráfico
  ggsave(
    filename = file.path(output_dir, paste0(gene, "_plot_prop_diplotipos.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message(paste("Procesado:", gene, "-> guardado en", output_dir))
}