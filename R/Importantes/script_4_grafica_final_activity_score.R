# =========================
# CYP2C9 + CYP2D6 Activity Score — COMPARISON
# =========================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

setwd("Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR")

clinpgx <- read_excel("frecuencias_poblacionales_activity_clinpgx.xlsx")
geno <- read_excel("salida_genostar_filtrado.xlsx")

cyps <- c("CYP2C9", "CYP2D6")

# =========================================================
# FUNCTION — OUR DATA
# =========================================================

procesar_cyp <- function(df, cyp){
  
  colname <- paste0(cyp, "_AS_1")
  
  if(!(colname %in% colnames(df))) return(NULL)
  
  df %>%
    mutate(AS = as.character(.data[[colname]])) %>%
    filter(
      !is.na(AS),
      !AS %in% c(
        "NA","N/A","","nan",
        "Invalid",
        "Diplotype not found",
        "Allele not found"
      )
    ) %>%
    mutate(
      AS_num = suppressWarnings(as.numeric(str_trim(AS)))
    ) %>%
    count(AS_num, .drop = FALSE) %>%
    mutate(
      freq = n / sum(n) * 100,
      Source = "Our Data",
      CYP = cyp
    )
}

geno_df <- bind_rows(lapply(cyps, function(x) procesar_cyp(geno, x)))

# =========================================================
# CLINPGX
# =========================================================

clin_df <- clinpgx %>%
  pivot_longer(
    cols = -Activity,
    names_to = "CYP",
    values_to = "freq"
  ) %>%
  rename(AS = Activity) %>%
  mutate(
    AS_clean = str_remove(AS, "≥"),
    AS_num = suppressWarnings(as.numeric(AS_clean)),
    Source = "ClinPGx"
  ) %>%
  filter(CYP %in% cyps)

# =========================================================
# COMBINE
# =========================================================

all_data <- bind_rows(geno_df, clin_df)

# =========================================================
# PHENOTYPE RULES
# =========================================================

all_data <- all_data %>%
  mutate(
    
    phenotype = case_when(
      
      CYP == "CYP2C9" & AS_num %in% c(0, 0.5) ~ "Poor",
      CYP == "CYP2C9" & AS_num %in% c(1, 1.5) ~ "Intermediate",
      CYP == "CYP2C9" & AS_num == 2 ~ "Normal",
      
      CYP == "CYP2D6" & AS_num == 0 ~ "Poor",
      CYP == "CYP2D6" & AS_num == 0.5 ~ "Intermediate",
      CYP == "CYP2D6" & AS_num == 0.25 ~ "Intermediate",
      CYP == "CYP2D6" & AS_num == 1 ~ "Intermediate",
      CYP == "CYP2D6" & AS_num >= 1.25 & AS_num <= 2.25 ~ "Normal",
      CYP == "CYP2D6" & AS_num >= 2.5 & AS_num <= 6 ~ "Ultrarapid",
      
      TRUE ~ "Indeterminate"
    ),
    
    Source = factor(Source, levels = c("ClinPGx", "Our Data")),
    CYP = factor(CYP, levels = cyps)
  )

# =========================================================
# COLORS
# =========================================================

fill_colors <- c(
  
  "Poor (ClinPGx)"          = alpha("#FFD6A5", 0.4),
  "Poor (Our Data)"         = "#C97A2B",
  
  "Intermediate (ClinPGx)"  = alpha("#FFF3B0", 0.4),
  "Intermediate (Our Data)" = "#C9B23E",
  
  "Normal (ClinPGx)"        = alpha("#B7E4C7", 0.4),
  "Normal (Our Data)"       = "#2E7D32",
  
  "Ultrarapid (ClinPGx)"    = alpha("#A7C7E7", 0.4),
  "Ultrarapid (Our Data)"   = "#1565C0",
  
  "Indeterminate (ClinPGx)" = alpha("#D3D3D3", 0.4),
  "Indeterminate (Our Data)"= "#666666"
)

# =========================================================
# LABELS (FIXED ORDER SAFE)
# =========================================================

all_data <- all_data %>%
  mutate(
    fill_group = paste0(phenotype, " (", as.character(Source), ")"),
    freq_label = ifelse(freq > 0.055, sprintf("%.2f%%", freq), ""),
    
    AS_order = AS_num,
    
    AS_plot = case_when(
      is.na(AS_num) ~ "Indeterminate\n(ClinPGx)",
      CYP == "CYP2D6" & AS_num %in% c(0.75, 2.25, 4) ~ paste0(AS_num, "\n(ClinPGx)"),
      TRUE ~ as.character(AS_num)
    )
  )

# FIX ORDER (NO DUPLICATES)
order_levels <- all_data %>%
  distinct(AS_plot, AS_order) %>%
  arrange(AS_order) %>%
  pull(AS_plot)

all_data$AS_plot <- factor(all_data$AS_plot, levels = order_levels)

# =========================================================
# BACKGROUND COLORS FOR ACTIVITY SCORE REGIONS
# =========================================================

bg_data <- data.frame(
  
  CYP = c(
    rep("CYP2C9", 4),
    rep("CYP2D6", 5)
  ),
  
  xmin = c(
    
    0.5,
    4.5,
    8.5,
    19.5,
    
    0.5,
    1.5,
    5.5,
    10.5,
    19.5
  ),
  
  xmax = c(
    
    3.5,
    7.5,
    9.5,
    20.5,
    
    1.5,
    5.5,
    10.5,
    19.5,
    20.5
  ),
  
  phenotype = c(
    "Poor",
    "Intermediate",
    "Normal",
    "Indeterminate",
    
    "Poor",
    "Intermediate",
    "Normal",
    "Ultrarapid",
    "Indeterminate"
  ),
  
  fill = c(
    "#E89B63",
    "#E6C15A",
    "#6FBF8F",
    "#9E9E9E",
    
    "#E89B63",
    "#E6C15A",
    "#6FBF8F",
    "#5C8FD6",
    "#9E9E9E"
  )
)

# =========================================================
# PLOT
# =========================================================

p <- ggplot(all_data, aes(x = AS_plot, y = freq, fill = fill_group)) +
  
  geom_rect(
    data = bg_data,
    inherit.aes = FALSE,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = bg_data$fill,
    alpha = 0.15
  ) +
  
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "grey30",
    size = 0.2
  ) +
  
  geom_text(
    aes(label = freq_label),
    position = position_dodge(width = 0.75),
    vjust = -0.3,
    size = 3
  ) +
  
  facet_wrap(~ CYP, ncol = 1, scales = "free_x") +
  
  scale_fill_manual(values = fill_colors, drop = FALSE) +
  
  labs(
    title = "Activity Score Distribution — CYP2C9 vs CYP2D6",
    x = "Activity Score",
    y = "Frequency (%)",
    fill = "Phenotype"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

# =========================================================
# SAVE
# =========================================================

ggsave(
  "plot_final_activity_score.png",
  p,
  width = 22,
  height = 12,
  dpi = 600
)