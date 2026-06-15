library(readxl)

## 1. Cargar Mvalues epi (Mvalues de las 4cg significativos del DMR de CYP3A5)
df_metilacion <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGsig_DMR_CYP3A5_Mvalues.xlsx")
df_metilacion <- as.data.frame(df_metilacion)      # convertir tibble → data.frame
rownames(df_metilacion) <- df_metilacion$CpG_ID    # ahora sí puedes asignar rownames
df_metilacion$CpG_ID <- NULL            # eliminar la columna


## 2. Trasponer para que sea mas facil 
df_metilacion_t <- as.data.frame(t(df_metilacion))

## 3. Calcular la media de los Mvalues de las 4cg para cda inidividuo 
df_metilacion_t$Mean_CpG <- rowMeans(df_metilacion_t[, 1:4], na.rm = TRUE)

## 4. Calcualr la mediana de las medias de todos los individuos 
mediana <- median(df_metilacion_t$Mean_CpG, na.rm = TRUE)

## 5.Clasificar los indivudios en metilacion alta y baja segun la mediana 
df_metilacion_t$Nivel_metilacion <- ifelse(df_metilacion_t$Mean_CpG > mediana, "alto", "bajo")
table(df_metilacion_t$Nivel_metilacion) # Comprobar que esta mitad y mitad 


## 6. Cargar transcriptomica gene_counts normalizados
df_trans <- read_excel("X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_epi_238_individuos.xlsx")
df_trans <- as.data.frame(df_trans)
df_trans$EnsemblID <- NULL #eLIMINAR COLUMNA DEL ENSEMBL ID 
df_trans <- df_trans[df_trans$CYP_name == "CYP3A5", ] #eLIMINAR TODAS LAS FILAS MENOS CYP3A5
rownames(df_trans) <- df_trans$CYP_name    # ahora sí puedes asignar rownames
df_trans$CYP_name <- NULL            # eliminar la columna


## 7. Trasponer para que sea mas facil 
df_trans_t <- as.data.frame(t(df_trans))

## 8. Ordenar individuos transcriptomica como epi 
df_trans_t <- df_trans_t[rownames(df_metilacion_t), , drop = FALSE]

## 9. Anadir la columna de gene_counts al data frame completo 
df_metilacion_t$CYP3A5_exp <- df_trans_t[, "CYP3A5"]

## 10.Boxplot de todas las CpGs 
library(ggplot2)

ggplot(df_metilacion_t, aes(x = Nivel_metilacion, y = CYP3A5_exp, fill = Nivel_metilacion)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("bajo" = "#56B4E9", "alto" = "#E69F00")) +
  labs(
    title = "CYP3A5 normalized counts",
    x = "Nivel de metilación",
    y = "Normalized counts"
  ) +
  theme_minimal(base_size = 14)


# Guardar el gráfico
ggsave(
  "X:/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/Boxplot_DMR_CYP3A5_all_CpGs.png",
  width = 6,
  height = 5,
  dpi = 300
)

