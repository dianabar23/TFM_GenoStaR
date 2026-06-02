################################################################################
# SETUP
################################################################################

# LOAD PATHS
source("2_Paths.R")

# LAMBDA FUNCTION
calculate_lambda <- function(pvalues){
  chisq <- qchisq(1-pvalues, 1)
  lambda <- median(chisq)/qchisq(0.5, 1)
  return(lambda)
}



################################################################################
# LOAD LIBRARIES & MAKE FOLDERS TO EXPORT
################################################################################

suppressMessages({
  # LOAD LIBRARIES
  library(readxl)
  library(qqman)
  library(ggplot2)
  library(parallel)
  library(doParallel)
  library(writexl)
  library(car)
  library(sva)
  library(limma)
})


# MAKE FOLDERS
if (!dir.exists("output/Figures/Sex_adjusted")) dir.create("output/Figures/Sex_adjusted")
if (!dir.exists("output/Tables/Sex_adjusted")) dir.create("output/Tables/Sex_adjusted")
if (!dir.exists("output/RData/Sex_adjusted")) dir.create("output/RData/Sex_adjusted")

################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# LOAD SAMPLESHEET & EUROPEAN INFO
samplesheet <- as.data.frame(read_excel("output/Tables/samplesheet_after_QC.xlsx"))

# LOAD M VALUES & BETA VALUES
load("output/RData/Mvalues_final_autosomes.RData")
load("output/RData/beta_final_autosomes.RData")

# CTRL PROBES
load("output/RData/Positive_ctrlprobe_intensities.RData") #ctrl probes

# CELL TYPE PROPORTIONS
load("output/RData/Cellproportions.RData") 

# SNPs FOR ANCESTRY PCs
load("output/RData/comb_SNPs.RData")

# GENOTYPING DATA
if (available_genotyping_data == "yes"){
  load(genotyping_path)
}

################################################################################
# GENERAL
################################################################################

# SELECT SAMPLES
if (multiple_samples == "yes"){
  samples <- samplesheet[which(samplesheet$Baseline_sample == "Y"),]
} else {
  samples <- samplesheet
}
sample_names <- intersect(samples$Basename, colnames(Mvalues))
samples <- samples[which(samples$Basename %in% sample_names), ]
print(paste0("Nr of individuals included: ", length(samples$Basename)))
write_xlsx(samples, path = "output/Tables/Sex_adjusted/Final_samplesheet_CaseCtrl.xlsx")

Mvalues <- Mvalues[, colnames(Mvalues) %in% samples$Basename]
beta <- beta[, colnames(beta) %in% samples$Basename]

# VARIABLES
Age <- scale(as.numeric(samples$Age))[, 1]
smoking <- scale(Mvalues["cg05575921", ])[, 1]
CaseCtrl <- factor(samples$Case_Control, levels = c("Control", "Case"))
Sex <- factor(samples$Sex)

# CELL TYPE PCs
if (array.type == "EPICv2"){
  cellcounts_df <- as.data.frame(cellcounts)
} else {
  cellcounts_df = cellcounts[[1]]
}

cellcounts_df = cellcounts_df[rownames(cellcounts_df) %in% sample_names, ]

cell_PCs = prcomp(cellcounts_df)
cell_PCs = as.data.frame(cell_PCs$x)
colnames(cell_PCs) <- paste0("cell_", colnames(cell_PCs))
cell_PCs <- sapply(cell_PCs, scale)



# CTRL PROBE PCs
ctrl <- ctrl[, colnames(ctrl) %in% sample_names]
pca <- prcomp(na.omit(t(ctrl))) # run pca
ctrlprobe_PCAscores <- as.data.frame(pca$x) #extract PCA scores
colnames(ctrlprobe_PCAscores) <- paste0("Ctrl_", colnames(ctrlprobe_PCAscores))
ctrlprobe_PCAscores <- sapply(ctrlprobe_PCAscores, scale)

#### ANCESTRY PCs ###
# PCA
if (available_genotyping_data == "yes"){
  
  if (!sample_names %in% colnames(genotype_matrix)){
    stop("PIPELINE STOPPED. Not all sample names are available in colnames(genotype_matrix). Please either ensure overlapping sample names (Basename in samplesheet) or use available_genotyping_data == no in 2_Paths.R")
  }
  
  genotype_matrix_red <- genotype_matrix[, colnames(genotype_matrix) %in% sample_names]
  pc <- prcomp(t(genotype_matrix_red))
}

if (available_genotyping_data == "no"){
  comb_SNPs_red <- comb_SNPs[, colnames(comb_SNPs) %in% sample_names]
  pc <- prcomp(t(comb_SNPs_red))
}

anc_PCs <- as.data.frame(pc$x)
colnames(anc_PCs) <- paste0("anc_", colnames(anc_PCs))
anc_PCs <- sapply(anc_PCs, scale)

# MAKE DF WITH VARIABLES FOR MODEL
variables_df <- data.frame(CaseCtrl = CaseCtrl, Sex = Sex, Age = Age, smoking = smoking)
variables_df <- cbind(variables_df, cell_PCs[,1:5], ctrlprobe_PCAscores[,1:15], anc_PCs[,1:15]) #he subido los Anc_PCs a 15
#variables_df <- variables_df[, !colnames(variables_df) %in% c("Ctrl_PC3", "Ctrl_PC4")] #he quitado los Ctrl_PCs 4
variables_df <- variables_df[, !colnames(variables_df) %in% c("Ctrl_PC4")] #se quita solo el Ctrl_Pc4

# TEST MULTI-COLINEARITY
temp_variables_df <- variables_df
temp_variables_df$Mvalues <- Mvalues[1, ]

if (nrow(samples) >= 35){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + #he quitado ctrl_PC 4
                                    Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                                    anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 + 
                                    anc_PC7 + anc_PC8 + anc_PC9 + anc_PC10 + anc_PC11 + anc_PC12 + anc_PC13 + anc_PC14 + anc_PC15, data = temp_variables_df))) #he añadido hasta el 15 Ancestry PC
  
}

if (nrow(samples) >= 25 & nrow(samples) < 35){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
                                    Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5, data = temp_variables_df)))
  
}

if (nrow(samples) < 25){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +
                                    anc_PC1 + anc_PC2, data = temp_variables_df)))
  
}

vif_obj$Variable <- rownames(vif_obj)
colnames(vif_obj) <- c("vif", "variable")
write_xlsx(vif_obj, path = "output/Tables/Sex_adjusted/VIF_CaseCtrl_Sex_adjusted.xlsx")
rm(temp_variables_df)

high_vif <- vif_obj[5:34, ]
high_vif <- high_vif$variable[which(high_vif$vif > 10)] #ya no 

################################################################################
# RUN MODEL
################################################################################

#Vamos a probar el modelo con 2 y con 3 ctrl probes (porque son los que más porcentaje de la varianza explican)
ctrlPCs_available <- grep("^Ctrl_PC", colnames(variables_df), value = TRUE)
ctrl_2 <- ctrlPCs_available[1:2] #dos primeros ctrl_probes
ctrl_3 <- ctrlPCs_available[1:3] #tres primeros ctrl_probes
ctrl_4 <- ctrlPCs_available[1:4] #cuatro primeros ctrl_probes (son el 1, el 2, el 3 y el 5)

#Vamos a probar el modelo con diferentes Anc PCs
ancPCs_available <- grep("^anc_PC", colnames(variables_df), value = TRUE)
anc_2  <- ancPCs_available[1:2]
anc_5  <- ancPCs_available[1:5]
anc_10 <- ancPCs_available[1:10]
anc_15 <- ancPCs_available[1:15]

#Usamos solo 2 Cell_PCs porque son los que explican la mayor parte de la varianza
cell_PCs <- c("cell_PC1", "cell_PC2")


#El resto de variables:
base_vars <- c("CaseCtrl", "Sex", "Age", "smoking") #son el resto de variables que no son PCs 


# SELECT M VALUES
Mvalues <- Mvalues[row.names(Mvalues) != "cg05575921", ] # remove cpg that is used to adjust for smoking

# TABLE TO SAVE ALL LAMBDA VALUES (tabla vacía para ir metiendo los lambda)
lambda_table <- data.frame(Model = c("2Ctrl_2Anc", "3Ctrl_2Anc", "4Ctrl_2Anc", "2Ctrl_5Anc", "3Ctrl_5Anc", "4Ctrl_5Anc",
                                     "2Ctrl_10Anc", "3Ctrl_10Anc", "4Ctrl_10Anc", "2Ctrl_15Anc", "3Ctrl_15Anc", "4Ctrl_15Anc"), lambda = NA)


##################
#Modelo 2Ctrl_2Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_2, anc_2)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[1,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_2CtrlPCs_2AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 3Ctrl_2Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_3, anc_2)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[2,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_3CtrlPCs_2AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 4Ctrl_2Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_4, anc_2)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[3,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_4CtrlPCs_2AncPCs.xlsx") # Guardar en Excel


##################
#Modelo 2Ctrl_5Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_2, anc_5)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[4,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_2CtrlPCs_5AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 3Ctrl_5Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_3, anc_5)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[5,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_3CtrlPCs_5AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 4Ctrl_5Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_4, anc_5)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[6,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_4CtrlPCs_5AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 2Ctrl_10Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_2, anc_10)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[7,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_2CtrlPCs_10AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 3Ctrl_10Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_3, anc_10)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[8,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_3CtrlPCs_10AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 4Ctrl_10Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_4, anc_10)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[9,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_4CtrlPCs_10AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 2Ctrl_15Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_2, anc_15)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_15AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_15AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[10,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_2CtrlPCs_15AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 3Ctrl_15Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_3, anc_15)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_15AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_3CtrlPCs_15AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[11,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_3CtrlPCs_15AncPCs.xlsx") # Guardar en Excel

##################
#Modelo 3Ctrl_15Anc
##################
model_vars <- c(base_vars, cell_PCs, ctrl_4, anc_15)
model_variables <- variables_df[, model_vars]
model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)
DMPs <- topTable(fit, num=Inf, coef=2)

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_15AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_4CtrlPCs_15AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[12,2] <- calculate_lambda(DMPs$P.Value)

# EXPORT
DMPs$CpG_ID <- rownames(DMPs)
DMPs <- DMPs[, c("CpG_ID", colnames(DMPs)[colnames(DMPs) != "CpG_ID"])] # Reordenar para que CpG_ID sea la primera columna
write_xlsx(DMPs, path = "output/Tables/Sex_adjusted/DMPs_BD_4CtrlPCs_15AncPCs.xlsx") # Guardar en Excel


# EXPORT LAMBDA TABLE
write_xlsx(lambda_table, path = "output/Tables/Sex_adjusted/Lambda_CaseCtrl_Sex_adjusted_limma.xlsx")


rm(list = ls());gc()
