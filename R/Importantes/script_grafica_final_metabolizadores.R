###
# =========================
# Frecuencias CYP - gráfico comparativo
# =========================

# Librerías
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# =========================
# 1. CARGA DE DATOS
# =========================
setwd("Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR")
clinpgx <- read_excel("frecuencias_poblacionales_metabolizadores_clinpgx.xlsx")
geno <- read_excel("salida_genostar_filtrado.xlsx")

# =========================
# 2. DEFINICIÓN CYPs A USAR
# =========================

cyps <- c("CYP3A5","CYP2B6","CYP2C9","CYP2C19","CYP2D6")

# Orden fijo de metabolizadores
met_order <- c(
  "Ultrarapid",
  "Rapid",
  "Normal",
  "Intermediate",
  "Poor",
  "Indeterminate",
  "Diplotype not found"
)

# =========================
# 3. FUNCIÓN NORMALIZACIÓN GENOSTAR
# =========================

procesar_cyp <- function(df, cyp){
  
  colname <- paste0(cyp, "_Metabolizer_Status")
  
  if(!(colname %in% colnames(df))){
    return(NULL)
  }
  
  data <- df %>%
    filter(!is.na(.data[[colname]])) %>%
    mutate(status = .data[[colname]])
  
  # Normalización de nombres (GENOSTAR -> categorías comunes)
  data <- data %>%
    mutate(status_clean = case_when(
      str_detect(status, "Ultra") ~ "Ultrarapid",
      str_detect(status, "Rapid") & !str_detect(status, "Ultra") ~ "Rapid",
      str_detect(status, "Normal") ~ "Normal",
      str_detect(status, "Intermediate") ~ "Intermediate",
      str_detect(status, "Poor") ~ "Poor",
      str_detect(status, "Invalid") ~ "Diplotype not found",
      str_detect(status, "Diplotype not found") ~ "Diplotype not found",
      str_detect(status, "Indeterminate") ~ "Indeterminate",
      TRUE ~ "Diplotype not found"
    ))
  
  out <- data %>%
    count(status_clean) %>%
    mutate(freq = n / sum(n) * 100,
           CYP = cyp,
           Source = "Our data") %>%
    select(CYP, status_clean, freq, Source)
  
  return(out)
}

# =========================
# 4. PROCESAR GENOSTAR
# =========================

geno_list <- lapply(cyps, function(cyp) procesar_cyp(geno, cyp))
geno_df <- bind_rows(geno_list)

# =========================
# 5. PROCESAR CLINPGX
# =========================

clin_df <- clinpgx %>%
  pivot_longer(
    cols = -Metabolizador,   # 👈 EXCLUYE la columna de texto
    names_to = "CYP",
    values_to = "freq"
  ) %>%
  rename(status_clean = Metabolizador) %>%
  mutate(
    Source = "ClinPGx"
  )

# =========================
# 6. UNIR DATOS
# =========================

all_data <- bind_rows(
  geno_df %>%
    rename(freq = freq),
  
  clin_df %>%
    rename(freq = freq)
)

# =========================
# 7. GRÁFICO
# =========================

# Asegurar orden de factores
all_data$status_clean <- factor(all_data$status_clean, levels = met_order)
all_data$Source <- factor(all_data$Source, levels = c("ClinPGx", "Our data"))
all_data$CYP <- factor(all_data$CYP, levels = cyps)


# Crear posición numérica para agrupar barras
all_data <- all_data %>%
  mutate(
    CYP_num = as.numeric(factor(CYP, levels = cyps)),
    offset = ifelse(Source == "ClinPGx", -0.2, 0.2),
    x_pos = CYP_num + offset
  )

p <- ggplot(all_data, aes(x = x_pos, y = freq, fill = status_clean, group = status_clean)) +
  
  geom_bar(
    stat = "identity",
    width = 0.35,
    color = "white",
    size = 0.2
  ) +
  
  # Etiquetas dentro de las barras
  geom_text(
    aes(label = ifelse(freq > 3, paste0(round(freq,1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  
  # Etiquetas debajo de cada barra (ClinPGx / Our data)
  geom_text(
    data = all_data %>%
      distinct(CYP, Source, x_pos),
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
      "Indeterminate" = "#D3D3D3",
      "Diplotype not found" = "#F5B7B1"
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


# =========================
# 8. GUARDAR FIGURA
# =========================

ggsave(
  "plot_final_metabolizadores.png",
  plot = p,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  "frecuencias_CYP_comparativa.pdf",
  plot = p,
  width = 14,
  height = 8
)

