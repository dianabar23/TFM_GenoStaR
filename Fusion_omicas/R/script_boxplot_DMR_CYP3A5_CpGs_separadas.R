library(readxl)



#### Mvalues (CON MEDIANA) frente a gene counts
## 1. Cargar Mvalues epi (Mvalues de las 4cg significativos del DMR de CYP3A5)
df_metilacion <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGsig_DMR_CYP3A5_Mvalues.xlsx")
df_metilacion <- as.data.frame(df_metilacion)      # convertir tibble → data.frame
rownames(df_metilacion) <- df_metilacion$CpG_ID    # ahora sí puedes asignar rownames
df_metilacion$CpG_ID <- NULL            # eliminar la columna


## 2. Trasponer para que sea mas facil 
df_metilacion_t <- as.data.frame(t(df_metilacion))


## 3. Calcualr la mediana de cada CpG de todos los individuos 
mediana_cg02084114 <- median(df_metilacion_t$cg02084114, na.rm = TRUE)
mediana_cg04706801 <- median(df_metilacion_t$cg04706801, na.rm = TRUE)
mediana_cg12694063 <- median(df_metilacion_t$cg12694063, na.rm = TRUE)
mediana_cg13596235 <- median(df_metilacion_t$cg13596235, na.rm = TRUE)

## 4.Clasificar los indivudios en metilacion alta y baja segun la mediana de cada CpG 
df_metilacion_t$Nivel_met_cg02084114 <- ifelse(df_metilacion_t$cg02084114 > mediana_cg02084114, "high", "low")
table(df_metilacion_t$Nivel_met_cg02084114) # Comprobar que esta mitad y mitad 

df_metilacion_t$Nivel_met_cg04706801 <- ifelse(df_metilacion_t$cg04706801 > mediana_cg04706801, "high", "low")
table(df_metilacion_t$Nivel_met_cg04706801) # Comprobar que esta mitad y mitad 

df_metilacion_t$Nivel_met_cg12694063 <- ifelse(df_metilacion_t$cg12694063 > mediana_cg12694063, "high", "low")
table(df_metilacion_t$Nivel_met_cg12694063) # Comprobar que esta mitad y mitad 

df_metilacion_t$Nivel_met_cg13596235 <- ifelse(df_metilacion_t$cg13596235 > mediana_cg13596235, "high", "low")
table(df_metilacion_t$Nivel_met_cg13596235) # Comprobar que esta mitad y mitad 



## 5. Cargar transcriptomica gene_counts normalizados
df_trans <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_epi_238_individuos.xlsx")
df_trans <- as.data.frame(df_trans)
df_trans$EnsemblID <- NULL #eLIMINAR COLUMNA DEL ENSEMBL ID 
df_trans <- df_trans[df_trans$CYP_name == "CYP3A5", ] #eLIMINAR TODAS LAS FILAS MENOS CYP3A5
rownames(df_trans) <- df_trans$CYP_name    # ahora sí puedes asignar rownames
df_trans$CYP_name <- NULL            # eliminar la columna


## 6. Trasponer para que sea mas facil 
df_trans_t <- as.data.frame(t(df_trans))

## 7. Ordenar individuos transcriptomica como epi 
df_trans_t <- df_trans_t[rownames(df_metilacion_t), , drop = FALSE]

## 8. Anadir la columna de gene_counts al data frame completo 
df_metilacion_t$CYP3A5_exp <- df_trans_t[, "CYP3A5"]


## 9. WILCOXON para significacion 
pvalor_wilcoxon_cg02084114 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg02084114, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg04706801 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg04706801, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg12694063 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg12694063, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg13596235 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg13596235, data = df_metilacion_t)$p.value

## 10. Boxplots con p-valores 
library(ggplot2)

p1 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg02084114, y = CYP3A5_exp, fill = Nivel_met_cg02084114)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg02084114",
    x = "Methylation level",
    y = "Normalized counts"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$CYP3A5_exp),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg02084114, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


p2 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg04706801, y = CYP3A5_exp, fill = Nivel_met_cg04706801)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg04706801",
    x = "Methylation level",
    y = "Normalized counts"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$CYP3A5_exp),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg04706801, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


p3 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg12694063, y = CYP3A5_exp, fill = Nivel_met_cg12694063)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg12694063",
    x = "Methylation level",
    y = "Normalized counts"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$CYP3A5_exp),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg12694063, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


p4 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg13596235, y = CYP3A5_exp, fill = Nivel_met_cg13596235)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg13596235",
    x = "Methylation level",
    y = "Normalized counts"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$CYP3A5_exp),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg13596235, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


library(patchwork)

figura_final <- (p1 | p2) / (p3 | p4)


# Guardar el gráfico
ggsave(
  "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/Boxplot_DMR_CYP3A5_CpGs_separadas_mediana.png",
  width = 15,
  height = 11,
  dpi = 300
)







#### Mvalues (CON MEDIANA) frente a Mvalues 
## 1. Cargar Mvalues epi (Mvalues de las 4cg significativos del DMR de CYP3A5)
df_metilacion <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGsig_DMR_CYP3A5_Mvalues.xlsx")
df_metilacion <- as.data.frame(df_metilacion)      # convertir tibble → data.frame
rownames(df_metilacion) <- df_metilacion$CpG_ID    # ahora sí puedes asignar rownames
df_metilacion$CpG_ID <- NULL            # eliminar la columna


## 2. Trasponer para que sea mas facil 
df_metilacion_t <- as.data.frame(t(df_metilacion))


## 3. Calcualr la mediana de cada CpG de todos los individuos 
mediana_cg02084114 <- median(df_metilacion_t$cg02084114, na.rm = TRUE)
mediana_cg04706801 <- median(df_metilacion_t$cg04706801, na.rm = TRUE)
mediana_cg12694063 <- median(df_metilacion_t$cg12694063, na.rm = TRUE)
mediana_cg13596235 <- median(df_metilacion_t$cg13596235, na.rm = TRUE)

## 4.Clasificar los indivudios en metilacion alta y baja segun la mediana de cada CpG 
df_metilacion_t$Nivel_met_cg02084114 <- ifelse(df_metilacion_t$cg02084114 > mediana_cg02084114, "alto", "bajo")
table(df_metilacion_t$Nivel_met_cg02084114) # Comprobar que esta mitad y mitad 

df_metilacion_t$Nivel_met_cg04706801 <- ifelse(df_metilacion_t$cg04706801 > mediana_cg04706801, "alto", "bajo")
table(df_metilacion_t$Nivel_met_cg04706801) # Comprobar que esta mitad y mitad 

df_metilacion_t$Nivel_met_cg12694063 <- ifelse(df_metilacion_t$cg12694063 > mediana_cg12694063, "alto", "bajo")
table(df_metilacion_t$Nivel_met_cg12694063) # Comprobar que esta mitad y mitad 

df_metilacion_t$Nivel_met_cg13596235 <- ifelse(df_metilacion_t$cg13596235 > mediana_cg13596235, "alto", "bajo")
table(df_metilacion_t$Nivel_met_cg13596235) # Comprobar que esta mitad y mitad 



## 9. WILCOXON para significacion 
pvalor_wilcoxon_cg02084114 <- wilcox.test(cg02084114 ~ Nivel_met_cg02084114, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg04706801 <- wilcox.test(cg04706801 ~ Nivel_met_cg04706801, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg12694063 <- wilcox.test(cg12694063 ~ Nivel_met_cg12694063, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg13596235 <- wilcox.test(cg13596235 ~ Nivel_met_cg13596235, data = df_metilacion_t)$p.value


## 10. Boxplots con p-valores 
library(ggplot2)

p1 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg02084114, y = cg02084114, fill = Nivel_met_cg02084114)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("bajo" = "#56B4E9", "alto" = "#E69F00")) +
  labs(
    title = "cg02084114",
    x = "Nivel de metilación",
    y = "M-value"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$cg02084114, na.rm = TRUE),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg02084114, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


p2 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg04706801, y = cg04706801, fill = Nivel_met_cg04706801)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("bajo" = "#56B4E9", "alto" = "#E69F00")) +
  labs(
    title = "cg04706801",
    x = "Nivel de metilación",
    y = "M-value"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$cg04706801, na.rm = TRUE),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg04706801, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



p3 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg12694063, y = cg12694063, fill = Nivel_met_cg12694063)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("bajo" = "#56B4E9", "alto" = "#E69F00")) +
  labs(
    title = "cg12694063",
    x = "Nivel de metilación",
    y = "M-value"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$cg12694063, na.rm = TRUE),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg12694063, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



p4 <- ggplot(df_metilacion_t, aes(x = Nivel_met_cg13596235, y = cg13596235, fill = Nivel_met_cg13596235)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("bajo" = "#56B4E9", "alto" = "#E69F00")) +
  labs(
    title = "cg13596235",
    x = "Nivel de metilación",
    y = "M-value"
  ) +
  annotate("text",
           x = 1.5,
           y = max(df_metilacion_t$cg13596235, na.rm = TRUE),
           label = paste0("p = ", signif(pvalor_wilcoxon_cg13596235, 3)),
           size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



library(patchwork)

figura_final <- (p1 | p2) / (p3 | p4)


# Guardar el gráfico
ggsave(
  "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/Boxplot_DMR_CYP3A5_CpGs_separadas_mediana_epi_epi.png",
  width = 15,
  height = 11,
  dpi = 300
)








###  Mvalues (CON DECILES) frente a gene countS
###  Mvalues (CON DECILES) frente a gene countS
## 1. Cargar Mvalues epi (Mvalues de las 4cg significativos del DMR de CYP3A5)
df_metilacion <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGsig_DMR_CYP3A5_Mvalues.xlsx")
df_metilacion <- as.data.frame(df_metilacion)      # convertir tibble → data.frame
rownames(df_metilacion) <- df_metilacion$CpG_ID    # ahora sí puedes asignar rownames
df_metilacion$CpG_ID <- NULL            # eliminar la columna


## 2. Trasponer para que sea mas facil 
df_metilacion_t <- as.data.frame(t(df_metilacion))


## 3. Calcualr la los deciles de cada individuo para cada CpG 
library(dplyr)
df_metilacion_t$decil_cg02084114 <- ntile(df_metilacion_t$cg02084114, 10)
df_metilacion_t$decil_cg04706801 <- ntile(df_metilacion_t$cg04706801, 10)
df_metilacion_t$decil_cg12694063 <- ntile(df_metilacion_t$cg12694063, 10)
df_metilacion_t$decil_cg13596235 <- ntile(df_metilacion_t$cg13596235, 10)



## 4.Clasificar los indivudios en metilacion alta y baja segun el decil (bajo para deciles 1,2,3 y 4 y alto para deciles 7,8,9 y 10)
df_metilacion_t$Nivel_met_cg02084114 <- case_when(
  df_metilacion_t$decil_cg02084114 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg02084114 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)

df_metilacion_t$Nivel_met_cg04706801 <- case_when(
  df_metilacion_t$decil_cg04706801 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg04706801 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)

df_metilacion_t$Nivel_met_cg12694063 <- case_when(
  df_metilacion_t$decil_cg12694063 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg12694063 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)

df_metilacion_t$Nivel_met_cg13596235 <- case_when(
  df_metilacion_t$decil_cg13596235 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg13596235 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)


## 5. Cargar transcriptomica gene_counts normalizados
df_trans <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_epi_238_individuos.xlsx")
df_trans <- as.data.frame(df_trans)
df_trans$EnsemblID <- NULL #eLIMINAR COLUMNA DEL ENSEMBL ID 
df_trans <- df_trans[df_trans$CYP_name == "CYP3A5", ] #eLIMINAR TODAS LAS FILAS MENOS CYP3A5
rownames(df_trans) <- df_trans$CYP_name    # ahora sí puedes asignar rownames
df_trans$CYP_name <- NULL            # eliminar la columna


## 6. Trasponer para que sea mas facil 
df_trans_t <- as.data.frame(t(df_trans))

## 7. Ordenar individuos transcriptomica como epi 
df_trans_t <- df_trans_t[rownames(df_metilacion_t), , drop = FALSE]

## 8. Anadir la columna de gene_counts al data frame completo 
df_metilacion_t$CYP3A5_exp <- df_trans_t[, "CYP3A5"]


## 9. WILCOXON para significacion 
pvalor_wilcoxon_cg02084114 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg02084114, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg04706801 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg04706801, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg12694063 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg12694063, data = df_metilacion_t)$p.value
pvalor_wilcoxon_cg13596235 <- wilcox.test(CYP3A5_exp ~ Nivel_met_cg13596235, data = df_metilacion_t)$p.value


library(dplyr)
library(ggplot2)

## 1. Construir dataframe largo
df_plot <- bind_rows(
  df_metilacion_t %>%
    select(CYP3A5_exp, Methylation_level = Nivel_met_cg04706801) %>%
    mutate(CpG = "cg04706801"),
  
  df_metilacion_t %>%
    select(CYP3A5_exp, Methylation_level = Nivel_met_cg12694063) %>%
    mutate(CpG = "cg12694063"),
  
  df_metilacion_t %>%
    select(CYP3A5_exp, Methylation_level = Nivel_met_cg13596235) %>%
    mutate(CpG = "cg13596235"),
  
  df_metilacion_t %>%
    select(CYP3A5_exp, Methylation_level = Nivel_met_cg02084114) %>%
    mutate(CpG = "cg02084114")
)

## 2. Filtrar NA (deciles 5–6)
df_plot <- df_plot %>% filter(!is.na(Methylation_level))

## 3. Tabla de p‑values
df_pvals <- data.frame(
  CpG = c("cg04706801", "cg12694063", "cg13596235", "cg02084114"),
  pval = c(
    pvalor_wilcoxon_cg04706801,
    pvalor_wilcoxon_cg12694063,
    pvalor_wilcoxon_cg13596235,
    pvalor_wilcoxon_cg02084114
  )
)

## 4. Gráfico único con colores pastel y p‑values pequeños
## Gráfico único con colores pastel, leyenda grande, CpGs más grandes y p‑values pequeños

p <- ggplot(df_plot, aes(x = CpG, y = CYP3A5_exp, fill = Methylation_level)) +
  geom_boxplot(alpha = 0.8, color = "black") +
  
  ## Colores pastel: low = rojo, high = verde
  scale_fill_manual(
    values = c(
      "low" = "#F4A6A6",   # rojo pastel
      "high" = "#A6E3A1"   # verde pastel
    )
  ) +
  
  ## Orden de las CpGs
  scale_x_discrete(limits = c("cg04706801", "cg12694063", "cg13596235", "cg02084114")) +
  
  ## p‑values pequeños
  geom_text(
    data = df_pvals,
    aes(
      x = CpG,
      y = max(df_plot$CYP3A5_exp) * 1.03,
      label = paste0("p = ", signif(pval, 3))
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  
  labs(
    title = "CYP3A5 Expression by CpG Methylation Level",
    x = NULL,
    y = "Normalized counts",
    fill = "Methylation level"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    ## CpGs más grandes y menos inclinadas
    axis.text.x = element_text(angle = 25, hjust = 1, size = 18),
    
    ## Leyenda más grande
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    
    ## Título más grande
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

print(p)


## 5. Guardar
ggsave(
  "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/Boxplot_DMR_CYP3A5_CpGs_separadas_deciles.png",
  p,
  width = 12,
  height = 8,
  dpi = 300
)






###  Mvalues (CON DECILES) frente a Mvalues 
## 1. Cargar Mvalues epi (Mvalues de las 4cg significativos del DMR de CYP3A5)
df_metilacion <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGsig_DMR_CYP3A5_Mvalues.xlsx")
df_metilacion <- as.data.frame(df_metilacion)      # convertir tibble → data.frame
rownames(df_metilacion) <- df_metilacion$CpG_ID    # ahora sí puedes asignar rownames
df_metilacion$CpG_ID <- NULL            # eliminar la columna


## 2. Trasponer para que sea mas facil 
df_metilacion_t <- as.data.frame(t(df_metilacion))


## 3. Calcualr la los deciles de cada individuo para cada CpG 
library(dplyr)
df_metilacion_t$decil_cg02084114 <- ntile(df_metilacion_t$cg02084114, 10)
df_metilacion_t$decil_cg04706801 <- ntile(df_metilacion_t$cg04706801, 10)
df_metilacion_t$decil_cg12694063 <- ntile(df_metilacion_t$cg12694063, 10)
df_metilacion_t$decil_cg13596235 <- ntile(df_metilacion_t$cg13596235, 10)



## 4.Clasificar los indivudios en metilacion alta y baja segun el decil (bajo para deciles 1,2,3 y 4 y alto para deciles 7,8,9 y 10)
df_metilacion_t$Nivel_met_cg02084114 <- case_when(
  df_metilacion_t$decil_cg02084114 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg02084114 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)

df_metilacion_t$Nivel_met_cg04706801 <- case_when(
  df_metilacion_t$decil_cg04706801 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg04706801 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)

df_metilacion_t$Nivel_met_cg12694063 <- case_when(
  df_metilacion_t$decil_cg12694063 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg12694063 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)

df_metilacion_t$Nivel_met_cg13596235 <- case_when(
  df_metilacion_t$decil_cg13596235 %in% 1:4 ~ "low",
  df_metilacion_t$decil_cg13596235 %in% 7:10 ~ "high",
  TRUE ~ NA_character_
)



## 10. Boxplots con p-valores 
library(ggplot2)

p1 <- ggplot(df_metilacion_t[!is.na(df_metilacion_t$Nivel_met_cg02084114), ],
             aes(x = Nivel_met_cg02084114, y = cg02084114, fill = Nivel_met_cg02084114)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg02084114",
    x = "Methylation level",
    y = "Mvalue"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


p2 <- ggplot(df_metilacion_t[!is.na(df_metilacion_t$Nivel_met_cg04706801), ],
         aes(x = Nivel_met_cg04706801, y = cg04706801, fill = Nivel_met_cg04706801))  +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg04706801",
    x = "Methylation level",
    y = "Mvalue"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



p3 <- ggplot(df_metilacion_t[!is.na(df_metilacion_t$Nivel_met_cg12694063), ],
             aes(x = Nivel_met_cg12694063, y = cg12694063, fill = Nivel_met_cg12694063)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg12694063",
    x = "Methylation level",
    y = "Mvalue"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



p4 <- ggplot(df_metilacion_t[!is.na(df_metilacion_t$Nivel_met_cg13596235), ],
             aes(x = Nivel_met_cg13596235, y = cg13596235, fill = Nivel_met_cg13596235)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("low" = "#56B4E9", "high" = "#E69F00")) +
  labs(
    title = "cg13596235",
    x = "Methylation level",
    y = "Mvalue"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



library(patchwork)

figura_final <- ((p1 | p2) / (p3 | p4)) +
  plot_annotation(
    title = "M-values of CpGs in CYP3A5 DMR",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )



# Guardar el gráfico
ggsave(
  "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/Boxplot_DMR_CYP3A5_CpGs_separadas_deciles_epi_epi.png",
  width = 15,
  height = 11,
  dpi = 300
)




