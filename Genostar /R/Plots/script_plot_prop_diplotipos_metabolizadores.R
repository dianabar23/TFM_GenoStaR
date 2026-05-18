
### CYP1A2 Y CYP2A4 no tienen _Diplotype_Phenotype_Table


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
library(readxl)
library(scales)

base_path <- "X:/Fobos/Proyecto_Diana/TFM_GenoStaR"

cyp_genes <- c("CYP3A5","CYP2C19","CYP2C9","CYP2B6")

for (gene in cyp_genes) {
  
  obj_name <- paste0(gene, "_resultados_final_diplotipos")
  
  if (!exists(obj_name)) {
    message(paste("No existe:", obj_name, "-> salto"))
    next
  }
  
  obj <- get(obj_name)
  
  output_dir <- file.path(base_path, gene)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---------------------------
  # DATOS DIPLOTIPOS
  # ---------------------------
  df <- obj[[1]] %>%
    select(LabID.V2, paste0(gene, "_diplotype")) %>%
    rename(Diplotype = 2)
  
  # ---------------------------
  # MAPA FENOTIPO DESDE EXCEL
  # ---------------------------
  pheno_file <- file.path(output_dir, paste0(gene, "_Diplotype_Phenotype_Table.xlsx"))
  
  pheno_map <- read_excel(pheno_file) %>%
    select(Diplotype = 1,
           Phenotype = 3)
  
  # ---------------------------
  # MERGE
  # ---------------------------
  df <- df %>%
    mutate(Diplotype = as.character(Diplotype)) %>%
    left_join(pheno_map, by = "Diplotype")
  
  # ---------------------------
  # FRECUENCIAS
  # ---------------------------
  freq <- df %>%
    count(Diplotype, Phenotype) %>%
    mutate(
      prop = n / sum(n),
      Phenotype = case_when(
        Diplotype == "NA" ~ "NA",
        TRUE ~ Phenotype
      ),
      Phenotype_group = case_when(
        Phenotype == "Poor" ~ "Poor",
        Phenotype %in% c("Intermediate", "Likely Intermediate") ~ "Intermediate",
        Phenotype == "Normal" ~ "Normal",
        Phenotype %in% c("Rapid", "Ultrarapid") ~ "Rapid/Ultrarapid",
        Phenotype == "Indeterminate" ~ "Indeterminate",
        Diplotype == "NA" ~ "NA",
        TRUE ~ "Other"
      )
    ) %>%
    arrange(desc(prop))
  
  # ---------------------------
  # GRAFICO
  # ---------------------------
  p <- ggplot(freq, aes(x = reorder(Diplotype, prop), y = prop, fill = Phenotype_group)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    geom_text(aes(label = percent(prop, accuracy = 0.1)),
              hjust = -0.1, size = 3) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(
      values = c(
        "NA" = "#D73027",
        "Poor" = "#542788",
        "Intermediate" = "#FEE08B",
        "Normal" = "#1A9850",
        "Rapid/Ultrarapid" = "#2C7BB6",
        "Indeterminate" = "#999999"
      )
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Distribución de diplotipos", gene),
      subtitle = "Coloreado por fenotipo farmacogenético",
      x = "Diplotipo",
      y = "Proporción",
      fill = "Fenotipo"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    expand_limits(y = max(freq$prop, na.rm = TRUE) * 1.15)
  
  # ---------------------------
  # GUARDAR
  # ---------------------------
  ggsave(
    filename = file.path(output_dir, paste0(gene, "_plot_prop_diplotipos_metabolizador.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message(paste("Procesado:", gene))
}