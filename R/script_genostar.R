#install.packages("devtools")
#devtools::install_github("GenoStaR-Genomics-Tools/GenoStaR", auth_token = NULL)

## Prueba con el documento de github 
library("GenoStaR")
setwd("C:/Users/CBM/Desktop/TFM_GenoStaR")
load("sample_CYP2D6_genotypes.rda")
genotypes <- sample_CYP2D6_genotypes
prueba <- all_geno_pheno(genotypes, c("CYP2D6"), phased = FALSE)


## GenoStaR con nuestros datos 
library(dplyr)
library("GenoStaR")
setwd("C:/Users/CBM/Desktop/TFM_GenoStaR")
load("matrix_geno_fixed_corregido.RData")
genotypes <- matrix_geno_fixed
salida <- all_geno_pheno(genotypes, 
                         c("CYP2D6","CYP1A2","CYP2B6","CYP2C9","CYP2C19",
                           "CYP3A4","CYP3A5"),
                         phased = FALSE, CYP1A2_name = "new")

## Edit Inés I 02_03_2026: proceso murió con error
#Error in rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep) : 
#  valor de 'times' no válido

## Prueba con nuestros datos por gen 
### 1. Se filtra por gen 

df_CYP1A2 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP1A2_rs")) 

df_CYP2B6 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2B6_rs")) 

df_CYP2C9 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2C9_rs")) 

df_CYP2C19 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2C19_rs")) 

df_CYP2D6 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2D6_rs")) 

df_CYP3A4 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A4_rs"))

df_CYP3A5 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A5_rs")) 

### 2. Se aplica el genostar a cada gen (para pocos individuos, con todos colapsa)
library("GenoStaR")
prueba_CYP3A5 <- all_geno_pheno(df_CYP3A5, c("CYP3A5"), phased = FALSE) #funciona

prueba_CYP3A4 <- all_geno_pheno(df_CYP3A4, c("CYP3A4"), phased = FALSE) #funciona

prueba_CYP1A2 <- all_geno_pheno(df_CYP1A2, c("CYP1A2"), phased = FALSE,  CYP1A2_name = "new") #funciona

prueba_CYP2D6 <- all_geno_pheno(df_CYP2D6, c("CYP2D6"), phased = FALSE)
prueba_CYP2D6 <- assign_diplotype(df_CYP2D6, "CYP2D6")

prueba_CYP2C19 <- all_geno_pheno(df_CYP2C19, c("CYP2C19"), phased = FALSE)
prueba_CYP2C19 <- assign_diplotype(df_CYP2C19, "CYP2C19")

prueba_CYP2C9 <- all_geno_pheno(df_CYP2C9, c("CYP2C9"), phased = FALSE)

prueba_CYP2B6 <- all_geno_pheno(df_CYP2B6, c("CYP2B6"), phased = FALSE)



#todos los cyp que funcionan hasta el momento juntos sobre un individuo 
salida <- all_geno_pheno(genotypes[22,], 
                         c("CYP1A2","CYP3A4","CYP3A5"),
                         phased = FALSE, CYP1A2_name = "new") #funciona




# CYP3A4
# 1. BUCLE OPTIMIZADO POR TIEMPO 
# este es el que se ha corrido, pasa al siguiente si tarda mas de 5 minutos 

library(dplyr)
library(R.utils)

df_CYP3A4 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A4_rs"))

n <- nrow(df_CYP3A4)

checkpoint_file <- "checkpoint_CYP3A4.rds"

# Si existe checkpoint, continuar desde ahí
if (file.exists(checkpoint_file)) {
  
  message("Cargando checkpoint...")
  resultados_lista <- readRDS(checkpoint_file)
  
  start_i <- which(sapply(resultados_lista, is.null))[1]
  
  if (is.na(start_i)) start_i <- n + 1
  
} else {
  
  resultados_lista <- vector("list", n)
  start_i <- 1
  
}

for (i in seq(from = start_i, to = n)) {
  
  message(paste0("Procesando muestra ", i, " de ", n))
  
  resultados_lista[[i]] <- tryCatch({
    
    withTimeout({
      
      resultado_temporal <- assign_diplotype(
        df_CYP3A4[i, ], 
        c("CYP3A4"), 
        phased = FALSE
      )
      
      resultado_temporal[[1]] %>%
        select(LabID.V2, CYP3A4_diplotype)
      
    }, timeout = 300)  # 5 minutos
    
  }, TimeoutException = function(ex) {
    
    message(paste("Timeout en muestra", i))
    
    data.frame(
      LabID.V2 = df_CYP3A4$LabID.V2[i],
      CYP3A4_diplotype = NA
    )
    
  }, error = function(e) {
    
    message(paste("Error en muestra", i, ":", e$message))
    
    data.frame(
      LabID.V2 = df_CYP3A4$LabID.V2[i],
      CYP3A4_diplotype = NA
    )
    
  })
  
  # liberar memoria
  gc()
  
  # guardar checkpoint cada 50 iteraciones
  if (i %% 50 == 0) {
    message("Guardando checkpoint...")
    saveRDS(resultados_lista, checkpoint_file)
  }
}
# resultado final
resultados_final <- bind_rows(resultados_lista)
# guardar resultado final
saveRDS(resultados_final, "CYP3A4_resultados_final.rds")



# 2. CREAR DATA FRAME CON LOS GENOTIPOS DE LOS NA 

df_CYP3A4 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A4_rs"))


df_CYP3A4_filtrado <- df_CYP3A4 %>%
  semi_join(
    CYP3A4_resultados_final %>%
      filter(is.na(CYP3A4_diplotype)),
    by = "LabID.V2"
  )

write.csv2(df_CYP3A4_filtrado, "df_CYP3A4_NA.csv", row.names = FALSE)



# CYP1A2 
# 1. BUCLE OPTIMIZADO POR TIEMPO 
# este es el que se ha corrido, pasa al siguiente si tarda mas de 5 minutos 

library(dplyr)
library(R.utils)

df_CYP1A2 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP1A2_rs"))

n <- nrow(df_CYP1A2)

checkpoint_file <- "checkpoint_CYP1A2.rds"

# Si existe checkpoint, continuar desde ahí
if (file.exists(checkpoint_file)) {
  
  message("Cargando checkpoint...")
  resultados_lista <- readRDS(checkpoint_file)
  
  start_i <- which(sapply(resultados_lista, is.null))[1]
  
  if (is.na(start_i)) start_i <- n + 1
  
} else {
  
  resultados_lista <- vector("list", n)
  start_i <- 1
  
}

for (i in seq(from = start_i, to = n)) {
  
  message(paste0("Procesando muestra ", i, " de ", n))
  
  resultados_lista[[i]] <- tryCatch({
    
    withTimeout({
      
      resultado_temporal <- assign_diplotype(
        df_CYP1A2[i, ], 
        c("CYP1A2"), 
        phased = FALSE, 
        CYP1A2_name = "new"
      )
      
      resultado_temporal[[1]] %>%
        select(LabID.V2, CYP1A2_diplotype)
      
    }, timeout = 300)  # 5 minutos
    
  }, TimeoutException = function(ex) {
    
    message(paste("Timeout en muestra", i))
    
    data.frame(
      LabID.V2 = df_CYP1A2$LabID.V2[i],
      CYP1A2_diplotype = NA
    )
    
  }, error = function(e) {
    
    message(paste("Error en muestra", i, ":", e$message))
    
    data.frame(
      LabID.V2 = df_CYP1A2$LabID.V2[i],
      CYP1A2_diplotype = NA
    )
    
  })
  
  # liberar memoria
  gc()
  
  # guardar checkpoint cada 50 iteraciones
  if (i %% 50 == 0) {
    message("Guardando checkpoint...")
    saveRDS(resultados_lista, checkpoint_file)
  }
}
# resultado final
resultados_final <- bind_rows(resultados_lista)
# guardar resultado final
saveRDS(resultados_final, "CYP1A2_resultados_final.rds")



# 2. CREAR DATA FRAME CON LOS GENOTIPOS DE LOS NA 

df_CYP1A2_filtrado <- df_CYP1A2 %>%
  semi_join(
    CYP1A2_resultados_final %>%
      filter(is.na(CYP1A2_diplotype)),
    by = "LabID.V2"
  )

write.csv2(df_CYP1A2_filtrado, "df_CYP1A2_NA.csv", row.names = FALSE)










# este no se ha corrido todavia, en vez de por tiempo es si llega a 10Gigas

library(dplyr)
library(pryr)

df_CYP1A2 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP1A2_rs"))
n <- nrow(df_CYP1A2)
checkpoint_file <- "checkpoint_CYP1A2.rds"
rm_limit <- 10 * 1024^3   # 10 GB en bytes

# cargar checkpoint si existe
if (file.exists(checkpoint_file)) {
  
  message("Cargando checkpoint...")
  resultados_lista <- readRDS(checkpoint_file)
  
  start_i <- which(sapply(resultados_lista, is.null))[1]
  
  if (is.na(start_i)) start_i <- n + 1
  
} else {
  
  resultados_lista <- vector("list", n)
  start_i <- 1
  
}

for (i in seq(from = start_i, to = n)) {
  
  message(paste0("Procesando muestra ", i, " de ", n))
  
  # comprobar uso de RAM
  current_ram <- pryr::mem_used()
  
  if (current_ram > ram_limit) {
    
    message("RAM > 10GB. Saltando muestra ", i)
    
    resultados_lista[[i]] <- data.frame(
      LabID.V2 = df_CYP1A2$LabID.V2[i],
      CYP1A2_diplotype = NA
    )
    
    gc()
    next
  }
  
  resultados_lista[[i]] <- tryCatch({
    
    resultado_temporal <- assign_diplotype(
      df_CYP1A2[i, ], 
      c("CYP1A2"), 
      phased = FALSE, 
      CYP1A2_name = "new"
    )
    
    resultado_temporal[[1]] %>%
      select(LabID.V2, CYP1A2_diplotype)
    
  }, error = function(e) {
    
    message(paste("Error en muestra", i, ":", e$message))
    
    data.frame(
      LabID.V2 = df_CYP1A2$LabID.V2[i],
      CYP1A2_diplotype = NA
    )
    
  })
  
  # limpiar objetos grandes
  rm(resultado_temporal)
  gc()
  
  # checkpoint cada 50
  if (i %% 50 == 0) {
    
    message("Guardando checkpoint...")
    saveRDS(resultados_lista, checkpoint_file)
    
  }
}

resultados_final <- bind_rows(resultados_lista)

saveRDS(resultados_final, "CYP1A2_resultados_final.rds")

