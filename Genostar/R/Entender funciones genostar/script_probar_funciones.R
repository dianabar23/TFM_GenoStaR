library(GenoStaR)
library(readxl)


library(GenoStaR)
library(readxl)
library(dplyr)
library(writexl)


PROBAR <- assign_diplotype_diana(df, c("CYP2C19"))
probar <- PROBAR[[1]]


#ruta <- "/Users/dianabarraso/Desktop/TFM_GenoStaR"
ruta <- "Y:/ctoma/Fobos/Proyecto_Diana/TFM_GenoStaR/Genostar"

# Leer archivos
df <- read_excel(file.path(ruta, "MATRIZ_FINAL_COMPLETA.xlsx"))
df <- as.data.frame(df)

datos_faltantes <- find_missing_data_diana(df)

### FUNCION: ALL GENO PHENO 
df_todo <- assign_diplotype_diana(df, c(
  "CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6",
  "CYP2C19",
  "CYP2D6"), phased = FALSE
  , CYP1A2_name = "new"
)


df_final <- df_todo[[1]]
df_final <- as.data.frame(df_final)

### Cargar df filtrada 
df <- read_xlsx("MATRIZ_FINAL_TODOS_LOS_GENES_GenoStaR.xlsx")
df[is.na(df)] <- " " # 2. Sustituyes los NA por un espacio en blanco " "

### Cargar df raw 
setwd("/Users/dianabarraso/Desktop/TFM_GenoStaR")
load("matrix_geno_fixed_espacios.RData")
df_raw <- matrix_geno_fixed



### FUNCION: find_missing_data()
missing_data <- find_missing_data(df)
missing_data_raw <- find_missing_data(df_raw)


### FUNCION: fill_empty_cells()
setwd("/Users/dianabarraso/Desktop/TFM_GenoStaR/CYP2C19")
load("CYP2C19_Allele_def.rda")
CYP2C19_allele_def <- CYP2C19_Allele_def
CYP2C19_allele_def <- CYP2C19_allele_def[-1, ]
CYP2C19_allele_def_COMPLETA <- fill_empty_cells(CYP2C19_allele_def)
## Se rellena 



### FUNCION: find_diplotype(data, genes, CYP1A2_name = )
diplotipos_raw <- find_diplotype(df_raw, c("CYP1A2","CYP3A4", "CYP3A5", "CYP2B6", "CYP2C19","CYP2C9"), CYP1A2_name = "new")
diplotipos <- find_diplotype(df, c("CYP1A2","CYP3A4", "CYP3A5", "CYP2B6", "CYP2C19","CYP2C9"), CYP1A2_name = "new")


diplotipos_df_raw <- diplotipos_raw[[1]]
diplotipos_df <- diplotipos[[1]]

no_match <- filter_no_match(diplotipos_df, "CYP1A2")




find_diplotype_diana_CYP2C19 <- function(df, gene, CYP1A2_name = "new") {
  
  gene <- toupper(gene)
  total_score <- 0
  
  # =========================
  # 1. Cargar tabla de alelos
  # =========================
  data <- switch(gene,
                 "CYP2C19" = CYP2C19_Allele_def,
                 stop("Gene not supported in this function")
  )
  
  missing_geno <- find_missing_data_diana(df)
  
  # quitar header row
  data <- data[-1, ]
  star_alleles <- data[[1]]
  data <- data[, -1, drop = FALSE]
  
  # =========================
  # 2. SNPs del paciente
  # =========================
  gene_snps <- grep(paste0("^", gene, "_"), colnames(df), value = TRUE)
  
  if (length(gene_snps) == 0) {
    stop("No SNP columns found for gene: ", gene)
  }
  
  # rs limpio
  snps_clean <- gsub(".*_(rs[0-9]+|X[0-9A-Z]+)", "\\1", gene_snps)
  
  # IMPORTANTÍSIMO: usar rs reales del panel
  valid_snps <- intersect(gsub("^rs", "", snps_clean),
                          gsub("^rs", "", colnames(data)))
  
  if (length(valid_snps) == 0) {
    stop("No overlapping SNPs between df and allele definition table")
  }
  
  # columnas reales del data
  data_mat <- data[, paste0("rs", valid_snps), drop = FALSE]
  
  # =========================
  # 3. filtrar filas válidas
  # =========================
  valid_rows <- apply(data_mat, 1, function(x) any(!is.na(x) & x != ""))
  
  data_mat <- data_mat[valid_rows, , drop = FALSE]
  star_alleles <- star_alleles[valid_rows]
  
  if (nrow(data_mat) == 0) {
    stop("No valid allele definitions after filtering SNP overlap")
  }
  
  # =========================
  # 4. preparar genotipo paciente
  # =========================
  df_gene <- df[, gene_snps, drop = FALSE]
  
  diplotypes <- rep(NA_character_, nrow(df))
  
  # =========================
  # 5. loop por paciente
  # =========================
  for (row_idx in seq_len(nrow(df_gene))) {
    
    genotypes <- lapply(df_gene[row_idx, ], function(x) strsplit(x, " ")[[1]])
    
    best_alleles <- c(NA, NA)
    best_scores <- c(-Inf, -Inf)
    
    # =========================
    # 6. matching contra cada allele
    # =========================
    for (j in seq_len(nrow(data_mat))) {
      
      score <- 0
      
      for (k in seq_along(valid_snps)) {
        
        patient_gt <- genotypes[[k]]
        ref_gt <- unlist(strsplit(data_mat[j, k], " "))
        
        # match simple
        if (any(patient_gt %in% ref_gt)) {
          score <- score + 1
        } else {
          score <- score - 2
        }
      }
      
      # asignar mejor match por haplotipo
      if (score > best_scores[1]) {
        best_scores[2] <- best_scores[1]
        best_alleles[2] <- best_alleles[1]
        
        best_scores[1] <- score
        best_alleles[1] <- star_alleles[j]
        
      } else if (score > best_scores[2]) {
        best_scores[2] <- score
        best_alleles[2] <- star_alleles[j]
      }
    }
    
    # =========================
    # 7. filtro de seguridad
    # =========================
    if (best_scores[1] < length(valid_snps) * 0.6) {
      diplotypes[row_idx] <- "No matching star alleles found"
    } else {
      diplotypes[row_idx] <- paste(best_alleles[1], best_alleles[2], sep = "/")
    }
  }
  
  df[[paste0(gene, "_diplotype")]] <- diplotypes
  
  return(list(df, missing_geno, total_score))
}

find_diplotype_diana_CYP2C19 <- function(df, gene, CYP1A2_name = "new") {
  
  gene <- toupper(gene)
  total_score <- 0
  
  # =========================
  # 1. Cargar Allele definition
  # =========================
  data <- switch(gene,
                 "CYP2C19" = CYP2C19_Allele_def,
                 stop("Gene not supported")
  )
  
  missing_geno <- find_missing_data_diana(df)
  
  # eliminar header interno
  data <- data[-1, ]
  star_alleles <- data[[1]]
  data <- data[, -1, drop = FALSE]
  
  # =========================
  # 2. SNPs del paciente
  # =========================
  gene_snps <- grep(paste0("^", gene, "_"), colnames(df), value = TRUE)
  
  if (length(gene_snps) == 0) {
    stop("No SNP columns found for gene: ", gene)
  }
  
  # extraer rs reales de forma segura
  snps_clean <- sub(".*_(rs[0-9]+|X.*)", "\\1", gene_snps)
  snps_clean <- gsub("^rs", "", snps_clean)
  
  # SNPs válidos en Allele_def
  allele_snps <- gsub("^rs", "", colnames(data))
  
  valid_snps <- snps_clean[snps_clean %in% allele_snps]
  
  if (length(valid_snps) == 0) {
    stop("No overlapping SNPs between df and allele definition")
  }
  
  # matriz de referencia
  data_mat <- data[, paste0("rs", valid_snps), drop = FALSE]
  
  # =========================
  # 3. filtrar alelos con datos reales
  # =========================
  valid_rows <- apply(data_mat, 1, function(x) any(!is.na(x) & x != ""))
  
  data_mat <- data_mat[valid_rows, , drop = FALSE]
  star_alleles <- star_alleles[valid_rows]
  
  if (nrow(data_mat) == 0) {
    stop("No valid allele definitions after filtering")
  }
  
  # =========================
  # 4. genotipos paciente
  # =========================
  df_gene <- df[, gene_snps, drop = FALSE]
  
  diplotypes <- rep(NA_character_, nrow(df))
  
  # =========================
  # 5. loop pacientes
  # =========================
  for (row_idx in seq_len(nrow(df_gene))) {
    
    best_alleles <- c(NA, NA)
    best_scores <- c(-Inf, -Inf)
    
    genotypes <- lapply(df_gene[row_idx, ], function(x) {
      if (is.na(x) || x == "") return(character(0))
      strsplit(x, " ")[[1]]
    })
    
    # =========================
    # 6. scoring contra alelos
    # =========================
    for (j in seq_len(nrow(data_mat))) {
      
      score <- 0
      coverage <- 0
      
      for (k in seq_along(valid_snps)) {
        
        ref <- data_mat[j, k]
        
        # penalizar alelos sin información
        if (is.na(ref) || ref == "" || ref == "NA") {
          score <- score - 2
          next
        }
        
        coverage <- coverage + 1
        
        patient_gt <- genotypes[[k]]
        ref_gt <- strsplit(ref, " ")[[1]]
        
        if (length(patient_gt) == 0) {
          score <- score - 2
        } else if (any(patient_gt %in% ref_gt)) {
          score <- score + 2
        } else {
          score <- score - 1
        }
      }
      
      # descartar alelos con poca información real
      if (coverage < length(valid_snps) * 0.6) {
        next
      }
      
      # guardar top 2
      if (score > best_scores[1]) {
        best_scores[2] <- best_scores[1]
        best_alleles[2] <- best_alleles[1]
        
        best_scores[1] <- score
        best_alleles[1] <- star_alleles[j]
        
      } else if (score > best_scores[2]) {
        best_scores[2] <- score
        best_alleles[2] <- star_alleles[j]
      }
    }
    
    # =========================
    # 7. decisión final
    # =========================
    if (best_scores[1] == -Inf) {
      diplotypes[row_idx] <- "No matching star alleles found"
    } else {
      diplotypes[row_idx] <- paste(best_alleles[1], best_alleles[2], sep = "/")
    }
  }
  
  df[[paste0(gene, "_diplotype")]] <- diplotypes
  
  return(list(df, missing_geno, total_score))
}

assign_diplotype_diana_CYP2C19 <- function(df, genes, phased = FALSE, CYP1A2_name = "new") {
  # Convert genes to uppercase to ensure consistency
  genes <- toupper(genes)
  
  # Identify SNP columns (containing "rs" or "CNV" in their names) to format genotypes, X is included for the special case of *18 in CYP1A2 which does not have rsID
  snp_columns <- grep("rs|CNV|X", names(df), value = TRUE)
  df[snp_columns] <- apply(df[snp_columns], 2, function(col) {
    sapply(col, format_genotype_diana)
  })
  
  # Initialize a list to store results for each gene
  all_gene_results <- vector("list", nrow(df))
  
  # Loop through all input genes
  for (gene in genes) {
    # Prepare columns for SNPs, get all relevant snp columns for the specified gene
    snp_cols <- grep(paste0(gene, "_(rs|CNV|X)"), names(df), value = TRUE)
    diplotype_col <- paste0(gene, "_diplotype")
    alternate_col <- paste0(gene, "_alternate_diplotype")
    
    # Initialize diplotype columns if not already in df
    if (!diplotype_col %in% names(df)) df[[diplotype_col]] <- NA
    if (!alternate_col %in% names(df)) df[[alternate_col]] <- vector("list", nrow(df))
    
    # Loop through each row in the dataframe
    for (row_idx in 1:nrow(df)) {
      # Check if all genotypes for this gene are NA, if yes assign NA and skip to next row
      if (all(is.na(df[row_idx, snp_cols, drop = FALSE]))) {
        df[row_idx, diplotype_col] <- NA
        df[[alternate_col]][[row_idx]] <- NA
        all_gene_results[[row_idx]] <- df[row_idx, ]
        next
      }
      # Identify the patient ID column dynamically
      patient_id_column <- names(df)[grep("ID", names(df), ignore.case = TRUE)][1]
      
      # Check if the patient ID column was found
      if (is.na(patient_id_column)) {
        stop("No patient ID column found. Ensure the dataframe contains a column with 'ID' in its name.")
      }
      
      # Include the patient ID column in the list of columns to select
      columns_to_select <- unique(c(patient_id_column, snp_cols))
      
      # Subset the dataframe
      snp_df <- df[row_idx, columns_to_select, drop = FALSE]
      patient_id <- snp_df[[patient_id_column]]
      snp_df <- snp_df[, !grepl("\\.1$", names(snp_df))]
      
      # Generate haplotype combinations or use snp_df directly
      if (phased) {
        haplotype_combinations <- snp_df
      } else {
        haplotype_combinations <- generate_phased_combinations_diana(snp_df, snp_cols)
      }
      
      # Process haplotypes
      match_results <- data.frame(Diplotype = character(), Score = numeric())
      for (haplotype_idx in 1:nrow(haplotype_combinations)) {
        temp_df <- snp_df
        if (!phased) {
          for (snp_idx in 1:length(snp_cols)) {
            haplotype <- haplotype_combinations[haplotype_idx, snp_idx]
            temp_df[[snp_cols[snp_idx]]] <- haplotype
          }
        }
        
        if (gene == "CYP2D6") {
          # Define the original expected CNV columns
          cnv_cols <- c("CYP2D6_CNV5UTR", "CYP2D6_CNVInt6", "CYP2D6_CNVx9")
          
          # Check for original CNV columns first
          available_cnv_cols <- cnv_cols[cnv_cols %in% colnames(df)]
          
          # If no original CNV columns are found, check for CYP2D6_CNV
          if (length(available_cnv_cols) == 0 && "CYP2D6_CNV" %in% colnames(df)) {
            available_cnv_cols <- "CYP2D6_CNV"
          }
          
          if (length(available_cnv_cols) == 0) {
            # No CNV columns available
            warning("No CNV columns found in the dataframe. Proceeding without CNVs.")
          } else {
            
            # Extract CNV data for the matching patient
            cnv_data <- df[df[[patient_id_column]] == patient_id, available_cnv_cols, drop = FALSE]
            
            #format cnv data to match expected input in code
            if (nrow(cnv_data) > 0) {
              # Correct single-copy values (e.g., "2" becomes "2 2")
              cnv_data[] <- lapply(cnv_data, function(col) {
                ifelse(nchar(col) == 1, paste(col, col), col)
              })
              
              # Rename `CYP2D6_CNV` to `CYP2D6_CNVx9` if it's the only available column
              if (length(available_cnv_cols) == 1 && available_cnv_cols == "CYP2D6_CNV" && 
                  !"CYP2D6_CNVx9" %in% colnames(cnv_data)) {
                colnames(cnv_data) <- "CYP2D6_CNVx9"
              }
              
              # Ensure the updated formatted cnv cols in `cnv_data` replaces the ones in `temp_df`
              if ("CYP2D6_CNVx9" %in% colnames(temp_df) && "CYP2D6_CNVx9" %in% colnames(cnv_data)) {
                temp_df$CYP2D6_CNVx9 <- NULL  # Remove the existing column
              }
              if ("CYP2D6_CNVInt6" %in% colnames(temp_df) && "CYP2D6_CNVInt6" %in% colnames(cnv_data)) {
                temp_df$CYP2D6_CNVInt6 <- NULL  # Remove the existing column
              }
              if ("CYP2D6_CNV5UTR" %in% colnames(temp_df) && "CYP2D6_CNV5UTR" %in% colnames(cnv_data)) {
                temp_df$CYP2D6_CNV5UTR <- NULL  # Remove the existing column
              }
              # Add to temp_df
              temp_df <- cbind(temp_df, cnv_data)
            } else {
              warning("No matching rows for the patient in the CNV data.")
            }
          }
        }
        
        #call Find diplotype function
        find_diplotype_result <- find_diplotype_diana_CYP2C19(temp_df, gene, CYP1A2_name)
        
        #store results
        result_df <- find_diplotype_result[[1]]
        missing_geno <- find_diplotype_result[[2]]
        total_score <- find_diplotype_result[[3]]
        result_status <- result_df[[diplotype_col]]
        if (any(result_status != "No matching star alleles found")) {
          match_results <- rbind(match_results, data.frame(Diplotype = result_status, Score = total_score))
        }
      }
      
      # Select best match
      if (nrow(match_results) > 0) {
        max_score <- max(match_results$Score)
        best_matches <- match_results[match_results$Score == max_score, "Diplotype"]
        unique_matches <- unique(sapply(best_matches, function(dip) {
          components <- strsplit(dip, "/")[[1]]
          paste(sort(components), collapse = "/")
        }))
        df[row_idx, diplotype_col] <- unique_matches[1]
        if (length(unique_matches) > 1) {
          df[[alternate_col]][[row_idx]] <- unique_matches[-1]
        }
      } else {
        df[row_idx, diplotype_col] <- NA
        df[[alternate_col]][[row_idx]] <- NA
      }
      
      # Store the updated row's result in the list
      all_gene_results[[row_idx]] <- df[row_idx, ]
    }
  }
  
  # Combine the list into a single dataframe
  final_results <- do.call(rbind, all_gene_results)
  
  # Remove alternate diplotype columns where all values are NA or NULL
  for (gene in genes) {
    alternate_col <- paste0(gene, "_alternate_diplotype")
    if (all(sapply(final_results[[alternate_col]], function(x) is.null(x) || is.na(x)))) {
      final_results[[alternate_col]] <- NULL
    }
  }
  
  # Convert list columns to character strings before returning
  final_results <- convert_list_to_string_diana(final_results)
  for(gene in genes){
    if(gene == "CYP2D6"){
      
      #adjust the potential alternate diplotypes 
      if("CYP2D6_alternate_diplotype" %in% colnames(final_results) || 
         ("CYP2D6_diplotype" %in% colnames(final_results) && 
          any(final_results$CYP2D6_diplotype %in% c("*10/*4", "*4/*10")))){
        
        final_results <- adjust_diplotype_diana(final_results)
        
        #check again if alternate diplotype column is all NA, if yes then remove it
        if (all(sapply(final_results[["CYP2D6_alternate_diplotype"]], function(x) is.null(x) || is.na(x)))) {
          final_results[["CYP2D6_alternate_diplotype"]] <- NULL
        }
      }
    }
  }
  
  #reorder diplotypes so smallest number is first
  # Identify diplotype columns (columns containing "diplotype" in the name)
  diplotype_cols <- grep("diplotype", colnames(final_results), value = TRUE)
  
  # Apply the reordering function to all diplotype columns
  #final_results[diplotype_cols] <- lapply(final_results[diplotype_cols], reorder_diplotype)
  final_results[diplotype_cols] <- lapply(final_results[diplotype_cols], function(x) {
    # Convert to character if the column is factor
    x <- as.character(x)
    # Apply the reorder_diplotype function to each element
    sapply(x, reorder_diplotype_diana)
  })
  
  return(list(final_results, missing_geno))
}

prueba <- assign_diplotype_diana_CYP2C19(df, c(
  "CYP2C19"), phased = FALSE
  , CYP1A2_name = "new"
)

prueba <- prueba[[1]]
