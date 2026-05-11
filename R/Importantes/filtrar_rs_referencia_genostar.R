library("GenoStaR")
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)


#Cargar el csv con todos los individuos 
load("matrix_geno_fixed_espacios.RData")
genotypes <- matrix_geno_fixed





### PARA CYP1A2: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP1A2_Allele_def_2.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP1A2_Allele_def_2)[-1]

#2. Coger solo los rs de este CYP 
df_CYP1A2 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP1A2_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP1A2 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP1A2_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP1A2_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )

#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta 
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP1A2)[-1] %>% 
  str_remove("CYP1A2_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP1A2_Allele_def_rs_excluidos <- CYP1A2_Allele_def_2 %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )


# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP1A2_Allele_def_rs_excluidos, "CYP1A2_Allele_def_rs_excluidos.xlsx")








### PARA CYP3A4: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP3A4_Allele_def.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP3A4_Allele_def)[-1] 

#2. Coger solo los rs de este CYP 
df_CYP3A4 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A4_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP3A4 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP3A4_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP3A4_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )


#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP3A4)[-1] %>% 
  str_remove("CYP3A4_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP3A4_Allele_def_rs_excluidos <- CYP3A4_Allele_def %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )

# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP3A4_Allele_def_rs_excluidos, "CYP3A4_Allele_def_rs_excluidos.xlsx")






### PARA CYP3A5: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP3A5_Allele_def.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP3A5_Allele_Def)[-1] 

#2. Coger solo los rs de este CYP 
df_CYP3A5 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP3A5_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP3A5 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP3A5_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP3A5_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )

#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta 
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP3A5)[-1] %>% 
  str_remove("CYP3A5_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP3A5_Allele_def_rs_excluidos <- CYP3A5_Allele_Def %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )


# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP3A5_Allele_def_rs_excluidos, "CYP3A5_Allele_def_rs_excluidos.xlsx")







### PARA CYP2C19: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP2C19_Allele_def.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP2C19_Allele_def)[-1]

#2. Coger solo los rs de este CYP 
df_CYP2C19 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2C19_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP2C19 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP2C19_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP2C19_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )

#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta 
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP2C19)[-1] %>% 
  str_remove("CYP2C19_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP2C19_Allele_def_rs_excluidos <- CYP2C19_Allele_def %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )


# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP2C19_Allele_def_rs_excluidos, "CYP2C19_Allele_def_rs_excluidos.xlsx")








### PARA CYP2C9: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP2C9_Allele_def.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP2C9_Allele_def)[-1]

#2. Coger solo los rs de este CYP 
df_CYP2C9 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2C9_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP2C9 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP2C9_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP2C9_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )


#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta 
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP2C9)[-1] %>% 
  str_remove("CYP2C9_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP2C9_Allele_def_rs_excluidos <- CYP2C9_Allele_def %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )


# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP2C9_Allele_def_rs_excluidos, "CYP2C9_Allele_def_rs_excluidos.xlsx")








### PARA CYP2B6: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP2B6_Allele_def.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP2B6_Allele_def)[-1]

#2. Coger solo los rs de este CYP 
df_CYP2B6 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2B6_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP2B6 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP2B6_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP2B6_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )

#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta 
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP2B6)[-1] %>% 
  str_remove("CYP2B6_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP2B6_Allele_def_rs_excluidos <- CYP2B6_Allele_def %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )


# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP2B6_Allele_def_rs_excluidos, "CYP2B6_Allele_def_rs_excluidos.xlsx")









### PARA CYP2D6: sobre TODOS los individuos, luego se hace el genostar
#1. Tabla de refencia de genostar (descargada de su github)
setwd("X:/Fobos/Proyecto_Diana/TFM_GenoStaR")
load("CYP2D6_Allele_def.rda")
# Extraemos los nombres de las columnas de referencia (ej: "rs123", "rs456")
# Asumimos que la primera columna es el nombre del alelo y las demás son los SNPs
cols_referencia <- colnames(CYP2D6_Allele_def)[-1]

#2. Coger solo los rs de este CYP 
df_CYP2D6 <- genotypes %>% 
  select(LabID.V2, starts_with("CYP2D6_rs"))

#3. Transformación y Filtrado Dinámico
df_resultado <- df_CYP2D6 %>%
  # Pasamos a formato largo (Tidy data)
  pivot_longer(
    cols = starts_with("CYP2D6_"), 
    names_to = "nombre_completo", 
    values_to = "genotipo"
  ) %>%
  # Creamos la columna limpia para comparar (quitando el prefijo)
  mutate(rs_limpio = sub("CYP2D6_", "", nombre_completo)) %>%
  # FILTRO: Solo nos quedamos con los rs que existen en tu archivo de referencia
  filter(rs_limpio %in% cols_referencia) %>%
  # (Opcional) Si quieres volver al formato original de una fila por paciente:
  pivot_wider(
    id_cols = 1, # Asumiendo que la col 1 es el ID del paciente
    names_from = nombre_completo, 
    values_from = genotipo
  )

#4. Calcular limitaciones: ver que *alelos no se tienen en cuenta 
#  Extraer los nombres de las columnas de los individuos 
# Quitamos el prefijo para que coincidan con los nombres en la tabla de definiciones
cols_presentes <- colnames(df_CYP2D6)[-1] %>% 
  str_remove("CYP2D6_")

#  Identificar y filtrar las columnas que NO están en los individuos
# Usamos !any_of para seleccionar las columnas cuyos nombres NO están en la lista
CYP2D6_Allele_def_rs_excluidos <- CYP2D6_Allele_def %>%
  select(
    1,                       # Mantenemos la primera columna (nombres de alelos/estrellas)
    !any_of(cols_presentes)  # Seleccionamos todas las que NO coincidan con los individuos
  )


# Guardar el dataframe resultante en un archivo .xlsx
write_xlsx(CYP2D6_Allele_def_rs_excluidos, "CYP2D6_Allele_def_rs_excluidos.xlsx")




