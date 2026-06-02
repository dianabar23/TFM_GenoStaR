library(readxl)
library(writexl)

# Leer archivos
gene_counts_trans <- read_excel(
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_transcriptomica_238_individuos.xlsx"
)

relacion_nombres <- read_excel(
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Relacion_nombres_muestras_epi_trancriptomica.xlsx"
)

# Crear diccionario Transcriptomica -> Epi
dic <- setNames(
  relacion_nombres$Epi,
  relacion_nombres$Transcriptomica
)

# Reemplazar nombres
nuevos_nombres <- sapply(colnames(gene_counts_trans), function(x){
  
  if(x %in% names(dic)){
    dic[[x]]
  } else {
    x
  }
  
})

# Asignar nuevos nombres
colnames(gene_counts_trans) <- nuevos_nombres

# Guardar
write_xlsx(
  gene_counts_trans,
  "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Fusion_omicas/Transcriptomica/CYP_gene_counts_normalizadas_nombres_epi_238_individuos.xlsx"
)
