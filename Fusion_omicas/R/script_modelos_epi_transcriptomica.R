### LOAD LIBRARIES
suppressMessages({
  library(readxl)
  library(limma)
  library(ggplot2)
  library(writexl)
  library(car)
})

##################
# Modelo Mvalues ~ gene_counts para cada CYP 
##################

### 1. CARGAR MVALUES (EPIGENÉTICA)
epi <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGs_por_CYP_ventana100kb_con_Mvalues.xlsx"

cyp_genes <- c(
  "CYP2D6",
  "CYP2B6",
  "CYP1A2",
  "CYP3A4",
  "CYP3A5",
  "CYP2C9",
  "CYP2C19"
)

Mvalues_CYPS <- lapply(cyp_genes, function(gene) {
  
  df <- read_excel(epi, sheet = gene)
  
  ## eliminar columnas de anotación
  df <- df[, -c(2:6)]
  
  df
})

names(Mvalues_CYPS) <- cyp_genes


### 2. CARGAR GENE COUNTS
gene_counts <- read_excel(
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_epi_238_individuos.xlsx"
)

gene_counts <- gene_counts[, -2]  # quitar EnsemblID


### 3. REORDENAR AMBAS VARIABLES PARA MISMO ORDEN DE MUESTRAS 

common_samples <- Reduce(intersect, list(
  colnames(Mvalues_CYPS[[1]])[-1],
  colnames(gene_counts)
))

## reordenar expresión
gene_counts <- gene_counts[, c("CYP_name", common_samples)]

## reordenar metilación
Mvalues_CYPS <- lapply(Mvalues_CYPS, function(df) {
  df[, c(names(df)[1], common_samples), drop = FALSE]
})


### ✔ CHECK CRÍTICO DE ALINEACIÓN (NUEVO)
stopifnot(all(common_samples == colnames(Mvalues_CYPS[[1]])[-1]))
stopifnot(all(common_samples == colnames(gene_counts)[-1]))


### 4. ESCALAR GENE COUNTS
gene_counts_scaled <- as.data.frame(
  scale(gene_counts[, -1])
)

gene_counts_scaled <- cbind(
  CYP_name = gene_counts$CYP_name,
  gene_counts_scaled
)


### 5. CREAR MATRICES (limma usa matrices)

df_gene_counts <- as.matrix(gene_counts_scaled[, -1])
rownames(df_gene_counts) <- gene_counts_scaled$CYP_name


Mvalues_CYPS <- lapply(Mvalues_CYPS, function(df) {
  
  cpg_ids <- df[[1]]
  
  mat <- as.matrix(df[, -1])
  mode(mat) <- "numeric"
  
  rownames(mat) <- cpg_ids
  
  mat
})


### 6. APLICAR MODELO SIMPLE SOBRE TODOS LOS CYPS 

### CYP2D6
CYP2D6_expr <- as.numeric(df_gene_counts["CYP2D6", ])

dim(Mvalues_CYPS[["CYP2D6"]])
length(CYP2D6_expr)

design_CYP2D6 <- model.matrix(~ CYP2D6_expr)
fit_CYP2D6 <- lmFit(Mvalues_CYPS[["CYP2D6"]], design_CYP2D6)
fit_CYP2D6 <- eBayes(fit_CYP2D6)

DMPs_CYP2D6 <- topTable(fit_CYP2D6, num = Inf, coef = "CYP2D6_expr")
DMPs_sig_CYP2D6 <- DMPs_CYP2D6[DMPs_CYP2D6$adj.P.Val < 0.05, ]


### CYP3A5
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])

dim(Mvalues_CYPS[["CYP3A5"]])
length(CYP3A5_expr)

design_CYP3A5 <- model.matrix(~ CYP3A5_expr)
fit_CYP3A5 <- lmFit(Mvalues_CYPS[["CYP3A5"]], design_CYP3A5)
fit_CYP3A5 <- eBayes(fit_CYP3A5)


DMPs_CYP3A5 <- topTable(fit_CYP3A5, num = Inf, coef = "CYP3A5_expr")

# Guardarlo para poder ver DMRs
DMPs_CYP3A5_export <- data.frame(
  CpG = rownames(DMPs_CYP3A5),
  DMPs_CYP3A5,
  row.names = NULL
)
write_xlsx(
  DMPs_CYP3A5_export,
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/DMRs_CYP3A5.xlsx"
)


DMPs_sig_CYP3A5 <- DMPs_CYP3A5[DMPs_CYP3A5$adj.P.Val < 0.05, ]

# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5 <- DMPs_CYP3A5[
  DMPs_CYP3A5$P.Value < (0.05 / nrow(DMPs_CYP3A5)), 
]


### CYP3A4
CYP3A4_expr <- as.numeric(df_gene_counts["CYP3A4", ])

dim(Mvalues_CYPS[["CYP3A4"]])
length(CYP3A4_expr)

design_CYP3A4 <- model.matrix(~ CYP3A4_expr)
fit_CYP3A4 <- lmFit(Mvalues_CYPS[["CYP3A4"]], design_CYP3A4)
fit_CYP3A4 <- eBayes(fit_CYP3A4)

DMPs_CYP3A4 <- topTable(fit_CYP3A4, num = Inf, coef = "CYP3A4_expr")
DMPs_sig_CYP3A4 <- DMPs_CYP3A4[DMPs_CYP3A4$adj.P.Val < 0.05, ]

# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A4 <- DMPs_CYP3A4[
  DMPs_CYP3A4$P.Value < (0.05 / nrow(DMPs_CYP3A4)), ]


### CYP1A2
CYP1A2_expr <- as.numeric(df_gene_counts["CYP1A2", ])

dim(Mvalues_CYPS[["CYP1A2"]])
length(CYP1A2_expr)

design_CYP1A2 <- model.matrix(~ CYP1A2_expr)
fit_CYP1A2 <- lmFit(Mvalues_CYPS[["CYP1A2"]], design_CYP1A2)
fit_CYP1A2 <- eBayes(fit_CYP1A2)

DMPs_CYP1A2 <- topTable(fit_CYP1A2, num = Inf, coef = "CYP1A2_expr")
DMPs_sig_CYP1A2 <- DMPs_CYP1A2[DMPs_CYP1A2$adj.P.Val < 0.05, ]


### CYP2C19
CYP2C19_expr <- as.numeric(df_gene_counts["CYP2C19", ])

dim(Mvalues_CYPS[["CYP2C19"]])
length(CYP2C19_expr)

design_CYP2C19 <- model.matrix(~ CYP2C19_expr)
fit_CYP2C19 <- lmFit(Mvalues_CYPS[["CYP2C19"]], design_CYP2C19)
fit_CYP2C19 <- eBayes(fit_CYP2C19)

DMPs_CYP2C19 <- topTable(fit_CYP2C19, num = Inf, coef = "CYP2C19_expr")
DMPs_sig_CYP2C19 <- DMPs_CYP2C19[DMPs_CYP2C19$adj.P.Val < 0.05, ]


### CYP2C9
CYP2C9_expr <- as.numeric(df_gene_counts["CYP2C9", ])

dim(Mvalues_CYPS[["CYP2C9"]])
length(CYP2C9_expr)

design_CYP2C9 <- model.matrix(~ CYP2C9_expr)
fit_CYP2C9 <- lmFit(Mvalues_CYPS[["CYP2C9"]], design_CYP2C9)
fit_CYP2C9 <- eBayes(fit_CYP2C9)

DMPs_CYP2C9 <- topTable(fit_CYP2C9, num = Inf, coef = "CYP2C9_expr")
DMPs_sig_CYP2C9 <- DMPs_CYP2C9[DMPs_CYP2C9$adj.P.Val < 0.05, ]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP2C9 <- DMPs_CYP2C9[
  DMPs_CYP2C9$P.Value < (0.05 / nrow(DMPs_CYP2C9)), ]


### CYP2B6
CYP2B6_expr <- as.numeric(df_gene_counts["CYP2B6", ])

dim(Mvalues_CYPS[["CYP2B6"]])
length(CYP2B6_expr)

design_CYP2B6 <- model.matrix(~ CYP2B6_expr)
fit_CYP2B6 <- lmFit(Mvalues_CYPS[["CYP2B6"]], design_CYP2B6)
fit_CYP2B6 <- eBayes(fit_CYP2B6)

DMPs_CYP2B6 <- topTable(fit_CYP2B6, num = Inf, coef = "CYP2B6_expr")
DMPs_sig_CYP2B6 <- DMPs_CYP2B6[DMPs_CYP2B6$adj.P.Val < 0.05, ]



##################
# Modelos con covariables 
##################

### 0. EPI Y GENE COUNTS YA ESTA DE ANTES (Pasos 1-5 de modelo simple)

### 1. CARGAR COVARIABLES 

setwd("Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Archivos_necesarios")

samplesheet <- as.data.frame(read_excel("samplesheet_after_QC.xlsx")) # LOAD SAMPLESHEET & EUROPEAN INFO (edad, sexo y caso/control)
load("Mvalues_final_autosomes.RData") # LOAD M VALUES (para la variable smoking)
load("Cellproportions.RData") # CELL TYPE PROPORTIONS (proporciones celulares)
load("Positive_ctrlprobe_intensities.RData") # CTRL PROBES (probes de control positivo)
load("comb_SNPs.RData") # SNPs FOR ANCESTRY PCs (CpGs en SNPs para ancestría)


### 2. ORDENAR SAMPLES (ALINEADO CON common_samples)

samples <- samplesheet[
  match(common_samples, samplesheet$Basename),
]

stopifnot(all(samples$Basename == common_samples))


### 3. AJUSTAR COVARIABLES

Age <- scale(as.numeric(samples$Age))[, 1]
Sex <- factor(samples$Sex)
CaseCtrl <- factor(samples$Case_Control, levels = c("Control", "Case"))


### SMOKING: depende del nivel de metilacion de la cg05575921, como es continua se escala
Mvalues_smoking <- Mvalues[, common_samples, drop = FALSE]
stopifnot("cg05575921" %in% rownames(Mvalues_smoking))
smoking <- scale(as.numeric(Mvalues_smoking["cg05575921", ]))[, 1]


### CELL TYPE: Hay que sacar primero los PC y luego escalar porque son continuos
cellcounts_df <- as.data.frame(cellcounts)
cellcounts_df <- cellcounts_df[common_samples, , drop = FALSE]
stopifnot(nrow(cellcounts_df) > 0)
cell_PCs = prcomp(cellcounts_df)
cell_PCs = as.data.frame(cell_PCs$x)
colnames(cell_PCs) <- paste0("cell_", seq_len(ncol(cell_PCs)))
cell_PCs <- as.data.frame(scale(cell_PCs))


### CTRL PROBES : Hay que sacar primero los PC y luego escalar porque son continuos
ctrl <- as.data.frame(ctrl)
ctrl <- ctrl[, common_samples, drop = FALSE]
stopifnot(ncol(ctrl) > 0)
ctrl_PCs <- prcomp(t(ctrl))
ctrl_PCs <- as.data.frame(ctrl_PCs$x)
colnames(ctrl_PCs) <- paste0("Ctrl_", seq_len(ncol(ctrl_PCs)))
ctrl_PCs <- as.data.frame(scale(ctrl_PCs))


### ANCESTRY: Hay que sacar primero los PC y luego escalar porque son continuos
comb_SNPs_red <- as.data.frame(comb_SNPs)
comb_SNPs_red <- comb_SNPs_red[, common_samples, drop = FALSE]
stopifnot(ncol(comb_SNPs_red) > 0)
anc_PCs <- prcomp(t(comb_SNPs_red))
anc_PCs <- as.data.frame(anc_PCs$x)
colnames(anc_PCs) <- paste0("anc_", seq_len(ncol(anc_PCs)))
anc_PCs <- as.data.frame(scale(anc_PCs))


### 4. DATAFRAME FINAL DE COVARIABLES (TU ESTRUCTURA ORIGINAL)

variables_df <- data.frame(CaseCtrl = CaseCtrl, Sex = Sex, Age = Age, smoking = smoking)

variables_df <- cbind(
  variables_df,
  cell_PCs[, 1:2],
  ctrl_PCs[, 1:2],
  anc_PCs[, 1:2]
)



### 5. APLICAR MODELOS ESPECIFICOS PARA CYP3A5 
##################
# MODELO COMPLETO: Mvalues ~ gene_counts + age + sex + smoking + Caso/Control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_MC <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + CaseCtrl + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2, data = variables_df)
fit_CYP3A5_MC <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_MC)
fit_CYP3A5_MC <- eBayes(fit_CYP3A5_MC)
DMPs_CYP3A5_MC <- topTable(fit_CYP3A5_MC, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_MC <- DMPs_CYP3A5_MC[DMPs_CYP3A5_MC$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_MC <- DMPs_CYP3A5_MC[DMPs_CYP3A5_MC$P.Value < (0.05 / nrow(DMPs_CYP3A5_MC)), ]


##################
# M1: Mvalues ~ gene_counts + age + sex 
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M1 <- model.matrix(~ CYP3A5_expr + Age + Sex, data = variables_df)
fit_CYP3A5_M1 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M1)
fit_CYP3A5_M1 <- eBayes(fit_CYP3A5_M1)
DMPs_CYP3A5_M1 <- topTable(fit_CYP3A5_M1, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M1 <- DMPs_CYP3A5_M1[DMPs_CYP3A5_M1$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M1 <- DMPs_CYP3A5_M1[DMPs_CYP3A5_M1$P.Value < (0.05 / nrow(DMPs_CYP3A5_M1)), ]


##################
# M2: Mvalues ~ gene_counts + age + sex + smoking 
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M2 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking, data = variables_df)
fit_CYP3A5_M2 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M2)
fit_CYP3A5_M2 <- eBayes(fit_CYP3A5_M2)
DMPs_CYP3A5_M2 <- topTable(fit_CYP3A5_M2, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M2 <- DMPs_CYP3A5_M2[DMPs_CYP3A5_M2$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M2 <- DMPs_CYP3A5_M2[DMPs_CYP3A5_M2$P.Value < (0.05 / nrow(DMPs_CYP3A5_M2)), ]


##################
# M3: Mvalues ~ gene_counts + age + sex + smoking + Caso/Control
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M3 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + CaseCtrl, data = variables_df)
fit_CYP3A5_M3 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M3)
fit_CYP3A5_M3 <- eBayes(fit_CYP3A5_M3)
DMPs_CYP3A5_M3 <- topTable(fit_CYP3A5_M3, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M3 <- DMPs_CYP3A5_M3[DMPs_CYP3A5_M3$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M3 <- DMPs_CYP3A5_M3[DMPs_CYP3A5_M3$P.Value < (0.05 / nrow(DMPs_CYP3A5_M3)), ]



##################
# MODELO 4: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M4 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2, data = variables_df)
fit_CYP3A5_M4 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M4)
fit_CYP3A5_M4 <- eBayes(fit_CYP3A5_M4)
DMPs_CYP3A5_M4 <- topTable(fit_CYP3A5_M4, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M4 <- DMPs_CYP3A5_M4[DMPs_CYP3A5_M4$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M4 <- DMPs_CYP3A5_M4[DMPs_CYP3A5_M4$P.Value < (0.05 / nrow(DMPs_CYP3A5_M4)), ]



##################
# MODELO 5: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M5 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + cell_1 + cell_2, data = variables_df)
fit_CYP3A5_M5 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M5)
fit_CYP3A5_M5 <- eBayes(fit_CYP3A5_M5)
DMPs_CYP3A5_M5 <- topTable(fit_CYP3A5_M5, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M5 <- DMPs_CYP3A5_M5[DMPs_CYP3A5_M5$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M5 <- DMPs_CYP3A5_M5[DMPs_CYP3A5_M5$P.Value < (0.05 / nrow(DMPs_CYP3A5_M5)), ]




##################
# MODELO 6: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M6 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + Ctrl_1 + Ctrl_2, data = variables_df)
fit_CYP3A5_M6 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M6)
fit_CYP3A5_M6 <- eBayes(fit_CYP3A5_M6)
DMPs_CYP3A5_M6 <- topTable(fit_CYP3A5_M6, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M6 <- DMPs_CYP3A5_M6[DMPs_CYP3A5_M6$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M6 <- DMPs_CYP3A5_M6[DMPs_CYP3A5_M6$P.Value < (0.05 / nrow(DMPs_CYP3A5_M6)), ]



##################
# MODELO 7: Mvalues ~ gene_counts + age + sex + smoking + anc_1 + anc_2
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M7 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + anc_1 + anc_2, data = variables_df)
fit_CYP3A5_M7 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M7)
fit_CYP3A5_M7 <- eBayes(fit_CYP3A5_M7)
DMPs_CYP3A5_M7 <- topTable(fit_CYP3A5_M7, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M7 <- DMPs_CYP3A5_M7[DMPs_CYP3A5_M7$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M7 <- DMPs_CYP3A5_M7[DMPs_CYP3A5_M7$P.Value < (0.05 / nrow(DMPs_CYP3A5_M7)), ]



##################
# MODELO 8: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M8 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2, data = variables_df)
fit_CYP3A5_M8 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M8)
fit_CYP3A5_M8 <- eBayes(fit_CYP3A5_M8)
DMPs_CYP3A5_M8 <- topTable(fit_CYP3A5_M8, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M8 <- DMPs_CYP3A5_M8[DMPs_CYP3A5_M8$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M8 <- DMPs_CYP3A5_M8[DMPs_CYP3A5_M8$P.Value < (0.05 / nrow(DMPs_CYP3A5_M8)), ]




##################
# MODELO 9: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + caso/control
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M9 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + CaseCtrl, data = variables_df)
fit_CYP3A5_M9 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M9)
fit_CYP3A5_M9 <- eBayes(fit_CYP3A5_M9)
DMPs_CYP3A5_M9 <- topTable(fit_CYP3A5_M9, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M9 <- DMPs_CYP3A5_M9[DMPs_CYP3A5_M9$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M9 <- DMPs_CYP3A5_M9[DMPs_CYP3A5_M9$P.Value < (0.05 / nrow(DMPs_CYP3A5_M9)), ]


##################
# MODELO 10: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + cell_1
##################
CYP3A5_expr <- as.numeric(df_gene_counts["CYP3A5", ])
design_CYP3A5_M10 <- model.matrix(~ CYP3A5_expr + Age + Sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + cell_1, data = variables_df)
fit_CYP3A5_M10 <- lmFit(Mvalues_CYPS[["CYP3A5"]],design_CYP3A5_M10)
fit_CYP3A5_M10 <- eBayes(fit_CYP3A5_M10)
DMPs_CYP3A5_M10 <- topTable(fit_CYP3A5_M10, num = Inf, coef = "CYP3A5_expr")
DMPs_sig_CYP3A5_M10 <- DMPs_CYP3A5_M10[DMPs_CYP3A5_M10$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A5_M10 <- DMPs_CYP3A5_M10[DMPs_CYP3A5_M10$P.Value < (0.05 / nrow(DMPs_CYP3A5_M10)), ]






### 6. APLICAR MODELOS ESPECIFICOS PARA CYP3A4
##################
# MODELO COMPLETO: Mvalues ~ gene_counts + age + sex + smoking + Caso/Control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2
##################
CYP3A4_expr <- as.numeric(df_gene_counts["CYP3A4", ])
design_CYP3A4_MC <- model.matrix(~ CYP3A4_expr + Age + Sex + smoking + CaseCtrl + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2, data = variables_df)
fit_CYP3A4_MC <- lmFit(Mvalues_CYPS[["CYP3A4"]],design_CYP3A4_MC)
fit_CYP3A4_MC <- eBayes(fit_CYP3A4_MC)
DMPs_CYP3A4_MC <- topTable(fit_CYP3A4_MC, num = Inf, coef = "CYP3A4_expr")
DMPs_sig_CYP3A4_MC <- DMPs_CYP3A4_MC[DMPs_CYP3A4_MC$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A4_MC <- DMPs_CYP3A4_MC[DMPs_CYP3A4_MC$P.Value < (0.05 / nrow(DMPs_CYP3A4_MC)), ]


##################
# M1: Mvalues ~ gene_counts + age + sex 
##################
CYP3A4_expr <- as.numeric(df_gene_counts["CYP3A4", ])
design_CYP3A4_M1 <- model.matrix(~ CYP3A5_expr + Age + Sex, data = variables_df)
fit_CYP3A4_M1 <- lmFit(Mvalues_CYPS[["CYP3A4"]],design_CYP3A4_M1)
fit_CYP3A4_M1 <- eBayes(fit_CYP3A4_M1)
DMPs_CYP3A4_M1 <- topTable(fit_CYP3A4_M1, num = Inf, coef = "CYP3A4_expr")
DMPs_sig_CYP3A4_M1 <- DMPs_CYP3A4_M1[DMPs_CYP3A4_M1$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A4_M1 <- DMPs_CYP3A4_M1[DMPs_CYP3A4_M1$P.Value < (0.05 / nrow(DMPs_CYP3A4_M1)), ]



##################
# MODELO 6: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 
##################
CYP3A4_expr <- as.numeric(df_gene_counts["CYP3A4", ])
design_CYP3A4_M6 <- model.matrix(~ CYP3A4_expr + Age + Sex + smoking + Ctrl_1 + Ctrl_2, data = variables_df)
fit_CYP3A4_M6 <- lmFit(Mvalues_CYPS[["CYP3A4"]],design_CYP3A4_M6)
fit_CYP3A4_M6 <- eBayes(fit_CYP3A4_M6)
DMPs_CYP3A4_M6 <- topTable(fit_CYP3A4_M6, num = Inf, coef = "CYP3A4_expr")
DMPs_sig_CYP3A4_M6 <- DMPs_CYP3A4_M6[DMPs_CYP3A4_M6$adj.P.Val < 0.05,]
# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_CYP3A4_M6 <- DMPs_CYP3A4_M6[DMPs_CYP3A4_M6$P.Value < (0.05 / nrow(DMPs_CYP3A4_M6)), ]





##################
# Modelo Mvalues ~ gene_counts para cada ZSCAN25 (los gene counts son del ZSCAN y las CpGs) 
##################
### 1. CARGAR MVALUES (EPIGENÉTICA)
epi <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/CpGs_por_CYP_ventana100kb_con_Mvalues.xlsx"

cyp_genes <- c("CYP3A5")

Mvalues_CYPS <- lapply(cyp_genes, function(gene) {
  
  df <- read_excel(epi, sheet = gene)
  
  ## eliminar columnas de anotación
  df <- df[, -c(2:6)]
  
  df
})

names(Mvalues_CYPS) <- cyp_genes




### 2. CARGAR GENE COUNTS
gene_counts <- read_excel(
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/ZSCAN25_gene_counts_normalizadas_nombres_epi_238_individuos.xlsx"
)

gene_counts <- gene_counts[, -2]  # quitar EnsemblID


### 3. REORDENAR AMBAS VARIABLES PARA MISMO ORDEN DE MUESTRAS 

common_samples <- Reduce(intersect, list(
  colnames(Mvalues_CYPS[[1]])[-1],
  colnames(gene_counts)
))

## reordenar expresión
gene_counts <- gene_counts[, c("CYP_name", common_samples)]

## reordenar metilación
Mvalues_CYPS <- lapply(Mvalues_CYPS, function(df) {
  df[, c(names(df)[1], common_samples), drop = FALSE]
})


### ✔ CHECK CRÍTICO DE ALINEACIÓN (NUEVO)
stopifnot(all(common_samples == colnames(Mvalues_CYPS[[1]])[-1]))
stopifnot(all(common_samples == colnames(gene_counts)[-1]))


### 4. ESCALAR GENE COUNTS 
# En este caso no se escala porque solo hay un gen 


### 5. CREAR MATRICES (limma usa matrices)

#gene counts sin escalado 
df_gene_counts <- as.matrix(gene_counts[, -1])
rownames(df_gene_counts) <- gene_counts$CYP_name

#epi
Mvalues_CYPS <- lapply(Mvalues_CYPS, function(df) {
  
  cpg_ids <- df[[1]]
  
  mat <- as.matrix(df[, -1])
  mode(mat) <- "numeric"
  
  rownames(mat) <- cpg_ids
  
  mat
})



### 6. APLICAR MODELO SIMPLE SOBRE ZSCAN25

ZSCAN25_expr <- as.numeric(df_gene_counts["ZSCAN25", ])

dim(Mvalues_CYPS[["CYP3A5"]])
length(ZSCAN25_expr)

design_ZSCAN25 <- model.matrix(~ ZSCAN25_expr)
fit_ZSCAN25 <- lmFit(Mvalues_CYPS[["CYP3A5"]], design_ZSCAN25)
fit_ZSCAN25 <- eBayes(fit_ZSCAN25)


DMPs_ZSCAN25 <- topTable(fit_ZSCAN25, num = Inf, coef = "ZSCAN25_expr")

# Guardarlo para poder ver DMRs
DMPs_ZSCAN25_export <- data.frame(
  CpG = rownames(DMPs_ZSCAN25),
  DMPs_ZSCAN25,
  row.names = NULL
)
write_xlsx(
  DMPs_ZSCAN25_export,
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Epi/DMPs_ZSCAN25.xlsx")



DMPs_sig_ZSCAN25 <- DMPs_ZSCAN25[DMPs_ZSCAN25$adj.P.Val < 0.05, ]

# ✔ Bonferroni dinámico (cambio mínimo)
DMPs_Bonferroni_ZSCAN25 <- DMPs_ZSCAN25[
  DMPs_ZSCAN25$P.Value < (0.05 / nrow(DMPs_ZSCAN25)), 
]

