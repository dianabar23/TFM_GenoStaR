
#Cargar CYP1A2_resultados_final_diplotipos
#Cargar CYP3A4_resultados_final_diplotipos
#Cargar CYP3A5_resultados_final_diplotipos
#Cargar CYP2C19_resultados_final_diplotipos
#Cargar CYP2C9_resultados_final_diplotipos
#Cargar CYP2B6_resultados_final_diplotipos
#Cargar CYP2D6_resultados_final_diplotipos




library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(writexl)

base_path <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

cyp_genes <- c("CYP1A2","CYP3A4","CYP3A5","CYP2C19","CYP2C9","CYP2B6","CYP2D6")

for (gene in cyp_genes) {
  
  obj_name <- paste0(gene, "_resultados_final_diplotipos")
  
  if (!exists(obj_name)) {
    message(paste("No existe:", obj_name))
    next
  }
  
  obj <- get(obj_name)
  
  output_dir <- file.path(base_path, gene)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---------------------------
  # DIPLOTIPOS
  # ---------------------------
  df <- obj[[1]] %>%
    select(LabID.V2, paste0(gene, "_diplotype")) %>%
    rename(Diplotype = 2)
  
  # ---------------------------
  # SEPARAR ALELOS (INCLUYE NA)
  # ---------------------------
  alleles <- df %>%
    mutate(Diplotype = as.character(Diplotype)) %>%
    separate_rows(Diplotype, sep = "/") %>%
    mutate(
      Allele = str_trim(Diplotype),
      Allele = ifelse(is.na(Allele), "NA", Allele)
    ) %>%
    select(Allele)
  
  # ---------------------------
  # FRECUENCIAS
  # ---------------------------
  freq_alleles <- alleles %>%
    count(Allele) %>%
    mutate(
      prop = n / sum(n),
      group = case_when(
        Allele == "NA" ~ "NA",
        Allele == "*1" ~ "Normal (*1)",
        TRUE ~ "Otros"
      )
    ) %>%
    arrange(desc(prop))
  
  # ---------------------------
  # GRAFICO
  # ---------------------------
  p <- ggplot(freq_alleles,
              aes(x = reorder(Allele, prop),
                  y = prop,
                  fill = group)) +
    
    geom_bar(stat = "identity") +
    coord_flip() +
    
    geom_text(aes(label = percent(prop, accuracy = 0.1)),
              hjust = -0.1, size = 3) +
    
    scale_y_continuous(labels = percent_format()) +
    
    scale_fill_manual(
      values = c(
        "NA" = "#D73027",
        "Normal (*1)" = "#1A9850",
        "Otros" = "#4575B4"
      )
    ) +
    
    theme_minimal(base_size = 12) +
    
    labs(
      title = paste("Frecuencia de alelos", gene, "(*star alleles)"),
      x = "Alelo",
      y = "Proporción",
      fill = ""
    ) +
    
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    
    expand_limits(y = max(freq_alleles$prop) * 1.15)
  
  # ---------------------------
  # GUARDAR PLOT
  # ---------------------------
  ggsave(
    filename = file.path(output_dir, paste0(gene, "_plot_prop_star_allele.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message(paste("Procesado:", gene))
}