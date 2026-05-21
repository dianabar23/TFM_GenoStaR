# =========================
# Frecuencias CYP - gráfico comparativo (CLEAN VERSION)
# =========================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

setwd("Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Genostar")

clinpgx <- read_excel("frecuencias_poblacionales_metabolizadores_clinpgx.xlsx")
geno <- read_excel("salida_genostar_filtrado.xlsx")

cyps <- c("CYP3A5","CYP2B6","CYP2C9","CYP2C19","CYP2D6")

met_order <- c(
  "Ultrarapid",
  "Rapid",
  "Normal",
  "Intermediate",
  "Poor",
  "Indeterminate"
)

# =========================================================
# FUNCIÓN GENOSTAR (CLEAN)
# =========================================================

procesar_cyp <- function(df, cyp){
  
  colname <- paste0(cyp, "_Metabolizer_Status")
  
  if(!(colname %in% colnames(df))){
    return(NULL)
  }
  
  data <- df %>%
    mutate(status = as.character(.data[[colname]])) %>%
    
    # 🔥 CLEAN DIRECTO: eliminar todo lo no válido
    filter(
      !is.na(status),
      !status %in% c("NA", "N/A", "", "nan", "Invalid", "Diplotype not found")
    ) %>%
    
    mutate(status_clean = case_when(
      str_detect(status, "Ultra") ~ "Ultrarapid",
      str_detect(status, "Rapid") & !str_detect(status, "Ultra") ~ "Rapid",
      str_detect(status, "Normal") ~ "Normal",
      str_detect(status, "Intermediate") ~ "Intermediate",
      str_detect(status, "Poor") ~ "Poor",
      str_detect(status, "Indeterminate") ~ "Indeterminate",
      TRUE ~ NA_character_
    )) %>%
    
    # 🔥 eliminar cualquier cosa no asignable
    filter(!is.na(status_clean))
  
  out <- data %>%
    count(status_clean) %>%
    mutate(
      freq = n / sum(n) * 100,
      CYP = cyp,
      Source = "Our data"
    )
  
  return(out)
}

# =========================================================
# GENOSTAR
# =========================================================

geno_list <- lapply(cyps, function(cyp) procesar_cyp(geno, cyp))
geno_df <- bind_rows(geno_list)

# =========================================================
# CLINPGX (sin cambios estructurales)
# =========================================================

clin_df <- clinpgx %>%
  pivot_longer(
    cols = -Metabolizador,
    names_to = "CYP",
    values_to = "freq"
  ) %>%
  rename(status_clean = Metabolizador) %>%
  mutate(Source = "ClinPGx")

# =========================================================
# UNIÓN
# =========================================================

all_data <- bind_rows(
  geno_df,
  clin_df
)

# =========================================================
# FACTORES
# =========================================================

all_data$status_clean <- factor(all_data$status_clean, levels = met_order)
all_data$Source <- factor(all_data$Source, levels = c("ClinPGx", "Our data"))
all_data$CYP <- factor(all_data$CYP, levels = cyps)

# =========================================================
# POSICIÓN EN X
# =========================================================

all_data <- all_data %>%
  mutate(
    CYP_num = as.numeric(factor(CYP, levels = cyps)),
    offset = ifelse(Source == "ClinPGx", -0.2, 0.2),
    x_pos = CYP_num + offset
  )

# =========================================================
# PLOT
# =========================================================

p <- ggplot(all_data, aes(x = x_pos, y = freq, fill = status_clean, group = status_clean)) +
  
  geom_bar(
    stat = "identity",
    width = 0.35,
    color = "grey30",
    size = 0.15
  ) +
  
  geom_text(
    aes(label = ifelse(freq > 0, sprintf("%.2f%%", freq), "")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  
  geom_text(
    data = all_data %>% distinct(CYP, Source, x_pos),
    aes(x = x_pos, y = -5, label = Source),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = unique(all_data$CYP_num),
    labels = cyps
  ) +
  
  scale_fill_manual(
    values = c(
      "Ultrarapid" = "#A7C7E7",
      "Rapid" = "#C6DBEF",
      "Normal" = "#B7E4C7",
      "Intermediate" = "#FFF3B0",
      "Poor" = "#FFD6A5",
      "Indeterminate" = "#D3D3D3"
    ),
    drop = FALSE
  ) +
  
  coord_cartesian(ylim = c(-5, 100)) +
  
  labs(
    title = "Frecuencias de metabolizadores de las enzimas CYP",
    x = "CYP",
    y = "Frecuencia (%)",
    fill = "Estado metabolizador"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0),
    legend.position = "right"
  )

# =========================================================
# SAVE
# =========================================================

ggsave(
  "plot_final_metabolizadores_2.png",
  plot = p,
  width = 14,
  height = 8,
  dpi = 600
)

