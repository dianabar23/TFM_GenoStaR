# Cargar librerias
library(tidyverse)
library(dplyr)

# Working directory 
setwd("C:/Users/CBM/Desktop/TFM_GenoStaR/")

#Los bfile se transforman con plink a PED y MAP
# Leer archivos ped y map (el normal y el adapatado)
map <- data.table::fread("660indivCYPS_SNPs_SNVs_05_03_26.map", header = FALSE, stringsAsFactors = FALSE)
ped <- data.table::fread("660indivCYPS_SNPs_SNVs_05_03_26.ped", header = FALSE, stringsAsFactors = FALSE)
map_formato <- data.table::fread("map_660indivCYPS_SNPs_SNVs_05_03_26_formato_geno.txt", header = FALSE, stringsAsFactors = FALSE) 
#este archivo sigue la misma estructura del map normal, lo unico es que a la segunda columna (rs) se le añade el CYP_ antes para que siga el formato para GenoStaR

## Sobre el map file 
colnames(map_formato) <- c("CHR", "SNP", "CM", "POS") # Darle nombre a las columna del map, para poder llamarlas luego
num_snps <- nrow(map) #Cada SNP tiene 2 columnas (alelo1, alelo2)

## Sobre el ped file 
geno <- ped[, -(1:6)] ## Separar la parte de genotipos del PED (6 primeras columnas)
# Comprobación de integridad
if (ncol(geno) != num_snps * 2) {
  stop("El número de columnas en PED no coincide con el número de SNPs en MAP.")
}

#Juntar los genotipos de dos en dos (ej:AA)
geno <- as.data.frame(geno)

matrix_geno <- as.data.frame(matrix(nrow = nrow(geno), ncol = num_snps))

for (i in seq_len(num_snps)) {
  matrix_geno[, i] <- paste0(geno[[2*i - 1]], geno[[2*i]])
}


#Poner el nombre las columnas como CYP_rs
colnames(matrix_geno) <- map_formato$SNP

#Añadir el labID 
ID <- ped[,2] # Crear la columna solo con el ID del paciente (corresponde a la segunda columna del ped file)
matrix_geno_final <- cbind(LabID = ID, matrix_geno)

# Cuando haya missing genotype que ponga una raya
matrix_geno_fixed <- as.data.frame(matrix_geno_final, stringsAsFactors = FALSE)
cols <- colnames(matrix_geno_fixed) != "LabID"
matrix_geno_fixed[, cols] <- lapply(matrix_geno_fixed[, cols], function(col) {
  col <- as.character(col)
  
  # quitar espacios para trabajar
  col_trim <- gsub(" ", "", col)
  
  # solo genotipos de 2 caracteres
  mask <- nchar(col_trim) == 2
  
  # en esos, cambiar 0 por -
  col_trim[mask] <- gsub("0", "-", col_trim[mask])
  
  # devolver la versión sin espacios (tipo AG, A-, -A, --)
  col_trim
})

#Guararlo como RData
save(matrix_geno_fixed, file = "matrix_geno_fixed_corregido.RData")



