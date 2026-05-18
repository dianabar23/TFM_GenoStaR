# =========================================================
# STAR ALLELE FREQUENCIES — CUSTOM LAYOUT (CORREGIDO)
# =========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readxl)
library(forcats)
library(tidytext)

base_path <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR"

cyp_genes <- c(
  "CYP1A2",
  "CYP3A4",
  "CYP3A5",
  "CYP2B6",
  "CYP2C19",
  "CYP2C9",
  "CYP2D6"
)

geno <- read_excel(
  file.path(base_path, "salida_genostar_filtrado.xlsx")
)

# =========================================================
# FUNCTION — PROCESS OUR DATA (FIXED)
# =========================================================

process_our_data <- function(df, gene){
  
  colname <- paste0(gene, "_diplotype")
  
  if(!colname %in% colnames(df)) return(NULL)
  
  raw <- df %>%
    mutate(
      diplotype = as.character(.data[[colname]]),
      
      diplotype = case_when(
        is.na(diplotype) |
          diplotype %in% c("NA", "N/A", "", "nan") ~
          NA_character_,
        TRUE ~ diplotype
      )
    )
  
  alleles <- raw %>%
    filter(!is.na(diplotype)) %>%
    separate_rows(diplotype, sep = "/") %>%
    mutate(
      allele = str_trim(diplotype)
    )
  
  counts <- alleles %>%
    count(allele)
  
  out <- counts %>%
    mutate(
      freq = n / sum(n) * 100,
      CYP = gene,
      Source = "Our data"
    )
  
  return(out)
}

# =========================================================
# CLINPGX FUNCTIONS (SIN CAMBIOS)
# =========================================================

read_clin_single <- function(gene){
  
  path <- file.path(
    base_path,
    "CYPs",
    gene,
    paste0(gene, "_frequency_table.xlsx")
  )
  
  if(!file.exists(path)) return(NULL)
  
  read_excel(path) %>%
    mutate(CYP = gene) %>%
    select(Allele, European100, CYP) %>%
    rename(
      allele = Allele,
      freq = European100
    ) %>%
    mutate(
      freq = as.numeric(
        str_replace_all(
          str_replace_all(as.character(freq), "%", ""),
          ",",
          "."
        )
      ),
      Source = "ClinPGx"
    )
}

clin_special <- read_excel(
  file.path(
    base_path,
    "frecuencias_poblacionales_allele_star_CYP1A2_CYP3A4_clinpgx.xlsx"
  )
) %>%
  pivot_longer(
    cols = -Allele,
    names_to = "CYP",
    values_to = "freq"
  ) %>%
  rename(allele = Allele) %>%
  mutate(
    freq = as.numeric(freq),
    Source = "ClinPGx"
  )

clin_single <- bind_rows(
  lapply(
    setdiff(cyp_genes, c("CYP1A2","CYP3A4")),
    read_clin_single
  )
)

clin_data <- bind_rows(
  clin_special,
  clin_single
)

# =========================================================
# OUR DATA
# =========================================================

our_alleles <- bind_rows(
  lapply(cyp_genes, function(x){
    process_our_data(geno, x)
  })
)

# =========================================================
# COMBINE DATA
# =========================================================

all_data <- bind_rows(
  our_alleles,
  clin_data
)

all_data <- all_data %>%
  filter(allele != "Diplotype not found")

# =========================================================
# KEEP ONLY RELEVANT ALLELES
# =========================================================

keep_alleles <- bind_rows(
  our_alleles %>% distinct(CYP, allele),
  clin_data %>% filter(freq > 1) %>% distinct(CYP, allele)
) %>%
  distinct()

all_data <- all_data %>%
  semi_join(keep_alleles, by = c("CYP", "allele"))

# =========================================================
# AESTHETICS (MODIFICADO)
# =========================================================

our_keys <- our_alleles %>%
  distinct(CYP, allele) %>%
  mutate(key = paste(CYP, allele))

all_data <- all_data %>%
  mutate(
    
    CYP = factor(
      CYP,
      levels = c(
        "CYP1A2","CYP3A4"," ","CYP3A5",
        "CYP2B6","CYP2C19","CYP2C9","CYP2D6")
    ),
    
    group = case_when(
      allele == "*1" ~ "Reference",
      TRUE ~ "Other"
    ),
    
    allele_source = paste0(allele, "    —    ", Source),
    freq_label = sprintf("%.2f%%", freq),
    
    fill_group = case_when(
      Source == "ClinPGx" & !(paste(CYP, allele) %in% paste(our_alleles$CYP, our_alleles$allele)) ~
        "Only ClinPGx (not in our data)",
      group == "Reference" & Source == "Our data" ~ "*1 (Our Data)",
      group == "Reference" & Source == "ClinPGx" ~ "*1 (ClinPGx)",
      group == "Other" & Source == "Our data" ~ "Other star alleles (Our Data)",
      TRUE ~ "Other star alleles (ClinPGx)"
    )
  )

# =========================================================
# ORDER (MODIFICADO PARA INDEPENDENCIA POR CYP)
# =========================================================

all_data <- all_data %>%
  group_by(CYP, allele) %>%
  mutate(
    # Frecuencia en "Our data" para este alelo/CYP concreto
    our_freq_val = max(freq[Source == "Our data"], na.rm = TRUE),
    # Si no existe en "Our data", le damos valor -1 para que vaya al final
    our_freq_val = ifelse(is.infinite(our_freq_val), -1, our_freq_val)
  ) %>%
  ungroup() %>%
  mutate(
    # reorder_within permite ordenar independientemente dentro de cada faceta (CYP)
    allele_ordered = reorder_within(allele, our_freq_val, CYP)
  )

# =========================================================
# PLOT (MODIFICADO)
# =========================================================

plot <- ggplot(
  all_data,
  aes(x = allele_ordered, y = freq, fill = fill_group)
) +
  
  geom_bar(stat = "identity", width = 0.7, position = position_dodge2(width = 0.8)) +
  
  geom_text(
    aes(label = freq_label),
    position = position_dodge2(width = 0.8),
    hjust = -0.1,
    size = 2.8
  ) +
  
  # scale_x_reordered limpia los sufijos internos de reorder_within
  scale_x_reordered() +
  
  coord_flip() +
  
  # scales = "free" es obligatorio para que reorder_within funcione correctamente
  facet_wrap(~ CYP, ncol = 3, scales = "free", drop = FALSE) +
  
  scale_fill_manual(
    values = c(
      "*1 (Our Data)"   = "#2E7D32",
      "*1 (ClinPGx)"  = "#66BB6A",
      "Other star alleles (Our Data)" = "#1565C0",
      "Other star alleles (ClinPGx)"= "#90CAF9",
      "Only ClinPGx (not in our data)" = "#BDBDBD"
    )
  ) +
  
  expand_limits(
    y = max(all_data$freq, na.rm = TRUE) * 1.25
  ) +
  
  labs(
    title = "Star allele frequencies across CYP enzymes",
    x = "",
    y = "Frequency (%)",
    fill = ""
  ) +
  
  theme_classic(base_size = 12) +
  
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 7.5),
    
    legend.position = c(0.9, 0.1),
    legend.justification = c(0.9, 0.1),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1.2, "cm"),
    
    legend.background = element_rect(
      fill = alpha("white", 0.6),
      colour = NA
    )
  )

# =========================================================
# SAVE
# =========================================================

ggsave(
  file.path(base_path, "plot_final_allele_star.png"),
  plot,
  width = 18,
  height = 14,
  dpi = 300
)
