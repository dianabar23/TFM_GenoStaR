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
library(scales)

# CYPs de interés
cyp_genes <- c("CYP3A5","CYP2C19","CYP2C9","CYP2B6")

# Lista vacía para unir todo
all_data <- data.frame()

for (gene in cyp_genes) {
  
  obj_name <- paste0(gene, "_resultados_final_diplotipos")
  
  if (!exists(obj_name)) {
    message(paste("No existe:", obj_name))
    next
  }
  
  obj <- get(obj_name)
  
  df <- obj[[1]] %>%
    select(LabID.V2, paste0(gene, "_diplotype")) %>%
    rename(Diplotype = 2)
  
  # Mapa fenotipo
  pheno_file <- file.path(
    "X:/Fobos/Proyecto_Diana/TFM_GenoStaR",
    gene,
    paste0(gene, "_Diplotype_Phenotype_Table.xlsx")
  )
  
  pheno_map <- readxl::read_excel(pheno_file) %>%
    select(Diplotype = 1, Phenotype = 3)
  
  df <- df %>%
    mutate(Diplotype = as.character(Diplotype)) %>%
    left_join(pheno_map, by = "Diplotype") %>%
    mutate(
      Phenotype_group = case_when(
        is.na(Diplotype) | Diplotype == "NA" ~ "NA",
        Phenotype == "Poor" ~ "Poor",
        Phenotype %in% c("Intermediate","Likely Intermediate") ~ "Intermediate",
        Phenotype == "Normal" ~ "Normal",
        Phenotype %in% c("Rapid","Ultrarapid") ~ "Rapid/Ultrarapid",
        Phenotype == "Indeterminate" ~ "Indeterminate",
        TRUE ~ "Indeterminate"
      ),
      CYP = gene
    )
  
  all_data <- bind_rows(all_data, df)
}

# ---------------------------
# PROPORCIONES
# ---------------------------
plot_data <- all_data %>%
  count(CYP, Phenotype_group) %>%
  group_by(CYP) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ---------------------------
# ORDEN
# ---------------------------
plot_data$Phenotype_group <- factor(
  plot_data$Phenotype_group,
  levels = c(
    "Rapid/Ultrarapid",
    "Intermediate",
    "Poor",
    "Normal",
    "Indeterminate",
    "NA"
  )
)


# ---------------------------
# GRAFICO
# ---------------------------
ggplot(plot_data, aes(x = CYP, y = prop, fill = Phenotype_group)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
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
    title = "Distribución de fenotipos metabolizadores por CYP",
    x = "CYP",
    y = "Proporción",
    fill = "Fenotipo"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )


# ---------------------------
# GUARDAR GRAFICO 
# ---------------------------
ggsave(
  filename = "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/CYPs_plot_prop_metabolizador.png",
  width = 10,
  height = 6,
  dpi = 300
)