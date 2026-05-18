# =========================================================
# CYP CALL RATE
# =========================================================

library(dplyr)
library(readxl)

# =========================================================
# PATH
# =========================================================

base_path <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR"

geno <- read_excel(
  file.path(base_path, "salida_genostar_filtrado.xlsx")
)

cyp_genes <- c(
  "CYP1A2",
  "CYP3A4",
  "CYP3A5",
  "CYP2B6",
  "CYP2C19",
  "CYP2C9",
  "CYP2D6"
)

# =========================================================
# FUNCTION
# =========================================================

get_call_rate <- function(df, gene){
  
  colname <- paste0(gene, "_diplotype")
  
  if(!colname %in% colnames(df)){
    return(NULL)
  }
  
  x <- df[[colname]]
  
  total <- length(x)
  
  missing <- sum(
    is.na(x) |
      x %in% c("NA", "N/A", "", "nan")
  )
  
  call_rate <- (total - missing) / total * 100
  
  cat(
    paste0(
      gene,
      " tiene un call rate de ",
      sprintf("%.2f", call_rate),
      "%\n"
    )
  )
}

# =========================================================
# RUN
# =========================================================

invisible(
  lapply(cyp_genes, function(g){
    get_call_rate(geno, g)
  })
)