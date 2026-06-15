# Proyecto Diana
## 19/02/26

### Objetivo
Obtener los bfile de los 660 individuos con las variantes SNPs y SNVs de los genes CYP (CYP1A2, CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP3A4 y CYP3A5)

### Pasos
#### 1. Sacar de los 5374 individuos, los 722 individuos con info de RNA 
Se eliminan 62 porque no estan en la tabla de los 5374 (ascendencia no europea u otras razones)

Resultado: 660indv.txt (familia (columna 1), ID individuo (columna 2))

#### 2. Sacar la info necesaria de los CYPs  
cromosoma, start y end de UCSC Genome Browser (a start y end se le pone una VENTANA de 5000)

Se hace una para el GRCh37/hg19 para las SNPs y se le llama CYPSconventana19SNP
Se hace una para el GRCh38/hg38 para las SNVs y se le llama CYPSconventana38SNV

id del gen de https://ensembl.org/

Resultado: CYPSconventana.txt


#### 3. PLINK para sacar los SNPs de los CYPs de los 660 individuos   
Los bfile de input: 5374_indiv_586000_SNPs_hg38_020725.hg19.ch.fl.bgn.fam tiene 5374 individuos y 5374_indiv_586000_SNPs_hg38_020725.hg19.ch.fl.bgn.bim tiene 9060687 SNPs

```bash
module load plink/1.90
plink --bfile bfile_og/5374_indiv_586000_SNPs_hg38_020725.hg19.ch.fl.bgn --keep 660indv.txt --extract 'range' CYPSconventana19SNP.txt --make-bed --out 660indivCYPS_SNPs_05_03_26
```

Los bfile de output: 660indivCYPS_SNPs_19_02_26.fam tiene 660 individuos y 660indivCYPS_SNPs_19_02_26.bim tiene 943 SNPs 

#### 4. PLINK para sacar los SNVs de los CYPs de los 660 individuos 
Se adapta el archivo de los 660indv.txt para que siga la estructura coorecta (ID individuo (columna 1), ID individuo (columna 2)) como 660_indiv_4_SNVs.list

Los bfile de input: 5668_indiv_110742_SNVs_hg38_091224.fam tiene 5668 individuos y 5668_indiv_110742_SNVs_hg38_091224.bim tiene 110742 SNVs

```bash
plink --bfile bfile_SNVs/5668_indiv_110742_SNVs_hg38_091224 --keep 660_indiv_4_SNVs.list --extract 'range' CYPSconventana38SNV.txt --make-bed --out 660indivCYPS_SNVs_05_03_26
```
Los bfile de output: 660indivCYPS_SNVs_05_03_26.fam tiene 660 individuos y 660indivCYPS_SNVs_05_03_26.bim tiene 342 SNVs 

#### 5. PLINK para fusionar los bfile con SNPs y SNVs 
Se adapta el archivo 660indivCYPS_SNPs_19_02_26.fam para que tenga el mismo formato que 660indivCYPS_SNVs_05_03_26.fam (1º y 2º columna con el ID del individuo) y se fusionan los bfile de ambos 
```bash
plink --bfile 660indivCYPS_SNPs_19_02_26 --bmerge 660indivCYPS_SNVs_05_03_26.bed 660indivCYPS_SNVs_05_03_26.bim 660indivCYPS_SNVs_05_03_26.fam --make-bed --out 660indivCYPS_SNPs_SNVs_05_03_26
```
Los bfile de output: 660indivCYPS_SNPs_SNVs_05_03_26.fam tiene los 660 individuos y 660indivCYPS_SNPs_SNVs_05_03_26.bim las 1260 variantes, se repiten 25 variantes

### Resultados 
Tras filtrar por individuos con info de RNA y familia quedan 660 individuos y tras filtrar por las variantes en los genes CYP quedan 1260 variantes

## 24/02/26
### Objetivo
Leer el artículo de GenoStaR y el Readme de gitbub para comprender como funciona, el input y output
### Resultados
Esquema (en el cuaderno) del artículo de GenoStaR y todo preperado para empezar a descargar GenoStaR en R 

## 26/02/26
### Objetivo
Formatear los bfile para convertirlos en el input para GenoStaR
### Pasos
#### Descargar GenoStaR en R studio 
Esta en el script de R script_genostar
#### Hacer la prueba con los datos de muestra del github 
Esta en el script de R script_genostar
#### Hacer el archivo input con el formato para genostar 
##### 1. Pasar los bfile a ped y map file (son mas sencillos de usar) con plink
```bash
 plink --bfile 660indivCYPS_SNPs_SNVs_05_03_26 --out 660indivCYPS_SNPs_SNVs_05_03_26 --recode
```
Te saca 660indivCYPS_SNPs_SNVs_05_03_26.map y 660indivCYPS_SNPs_SNVs_05_03_26.ped

##### 2. Conseg el map file adapatado juntando el excel el CYP al que pertence el rs y el rs (CYP_rs)

Con el excel concatenar_CYPS se formatea el map file y se guarda como map_660indivCYPS_SNPs_SNVs_05_03_26_formato_geno.txt

##### 3. Meter map, ped y map adaptado al el script de R formatear_input_genostar para generar la matrix_geno_fixed_05_03 (formato de input correcto para genostar)

#### Aplicar genostar sobre los datos 
Esta en el script de R script_genostar, se aplica sobre matrix_geno_fixed_05_03
### Resultados
La descarga ha ido bien, la prueba tambien, se ha formateado correctamente el input y se ha dejado corriendo genostar con nuestros datos  


## 03/03/26
### Objetivo
Optimizar genostar para nuestros datos 
### Pasos
Esta todo en el script de R script_genostar, se aplica sobre matrix_geno_fixed
#### Aplicar genostar sobre los datos  
No funciona, se queda pillado durante horas 
#### Aplicar genostar sobre los datos para cada CYP  
No funciona, se queda pillado durante horas 
#### Aplicar genostar sobre los datos para cada CYP de solo los primeros 5 individuos 
Funciona para CYP1A2 (new), CYP3A4 y CYP3A5. Se decide seguir con estos y ya se vera luego que pasa con los otros 
#### Aplicar genostar sobre el primer individuo para los tres CYPs que funcionan 
Sale bien 
### Resultados
Genostar no saca info para CYP2D6, CYP2B6, CYP2C19 y CYP2C9. Pero si que funciona para CYP1A2 (new), CYP3A4 y CYP3A5 por separado para pocos individuos. Para el primer individuo sale cuando se aplican todos juntos. Asi que se decide hacer un bucle que pase por todos los individuos para estos tres CYPs, para que no aumente el tiempo de CPU y se colapse como ocurría. 

## 05/03/26
### Objetivo
Optimizar genostar para nuestros datos 
### Pasos
Esta todo en el script de R script_genostar, se aplica sobre matrix_geno_fixed
#### Bucle Aplicar genostar sobre los datos  
Se hace y se intenta optimizar el bucle, pero no funciona porque necesita mucha RAM 
#### Repetir lo del día 19/02/26 y 26/02/26 para los SNV sobre GRCh38/hg38
Nos dimos cuenta de que las SNVs estaban hechas para GRCh37/hg19, mientras que los datos que nos habían dado eran de GRCh38/hg38, así que se volvió a hacer todo sobre GRCh38/hg38. (Se deja lo que se hizo corregido sobre esos días)
### Resultados
Se vuelve a crear todo pero bien, para el próximo día seguir optimizando GenoStaR y la RAM

## 10/03/26
### Objetivo
Optimizar genostar para nuestros datos 
### Pasos
Esta todo en el script de R script_genostar, se aplica sobre matrix_geno_fixed
#### Bucle Aplicar genostar sobre los datos  
Se prueban varios bucles, haciendo que no se genere un data frame por cada individuo (para que no ocupe tanto) y que genere un único dataframe con solo el individuo y su diplotipo, pero sigue dando el problema de la RAM 
### Resultados
La RAM sigue dando problemas aún optimizando lo máximo posible el bucle

## 12/03/26
### Objetivo
Optimizar genostar para nuestros datos 
### Pasos
Esta todo en el script de R script_genostar, se aplica sobre matrix_geno_fixed
#### Bucle Aplicar genostar sobre los datos  
Se prueban varios bucles, haciendo que no se genere un data frame por cada individuo (para que no ocupe tanto) y que genere un único dataframe con solo el individuo y su diplotipo, pero sigue dando el problema de la RAM 
### Resultados
La RAM sigue dando problemas aún optimizando lo máximo posible el bucle



## 17/03/26
### Objetivo
Optimizar genostar para nuestros datos 
### Pasos CYP1A2
Esta todo en el script de R script_genostar, se aplica sobre matrix_geno_fixed
#### Bucle Aplicar genostar sobre los datos  
Se optimiza el bucle para que si se tira más de 5 minutos en una muestra, pase al siguiente y ponga NA para ese individuo 
#### Se corre el bucle para CYP1A2: assign_diplotypes
Genera el objeto R CYP1A2_resultados_final. Hay info de diplotipo para 609 individuos y 51 no tiene diplotipo (NA)

Se sacan los genotipos de esos 51 individuos en df_CYP1A2_NA para poder trabajar con ellos 
#### Imputar los diplotipos de los individuos NA de CYP1A2  
Se hace a mano comparando el data frame con los genotipos de esos 51 individuos df_CYP1A2_NA y la tabla de internet CYP1A2_Haplotypes-PS216394-1454147960 que se filtra a CYP1A2_Haplotypes_FILTRADO (esta solo contiene los rs que coindicen para que sea más sencillo imputar)

### Pasos CYP3A4
#### Se corre el bucle de genostar para CYP3A4: assign_diplotypes
Se deja corriendo el bucle de script_genostar para CYP3A4
### Resultados
Para CYP1A2 ya se tienen los diplotipos de 609 individuos con genostar y 51 se intenta imputar. Se deja corriendo el bucle para CYP3A4


## 24/03/26
### Objetivo
Terminar de sacar los diplotipos de los 660 individuos para CYP1A2 y CYP3A4 
### Pasos CYP1A2 
#### Filtrar por los rs que usa genostar
Antes de seguir imputando a mano, se nos ocurre ver si el problema viene de que hay muchos rs en nuestros datos de inicio que Genostar no usa, para ello filtramos de los 51 individuos (df_CYP1A2_NA.csv) solo los rs que aparecern en la tabla de referencia que usa genostar, sacada de github genostar (CYP1A2_Allele_def_2.rda). Todo esta en el script filtrar_rs_referencia_genostar.
#### Correr el genostar de los NA: assign_diplotype 
Se ve que haciendo el filtrado se saca el diplotipo de los 51 pacientes restantes (CYP1A2_resultados_final_NA.rds), por loq eu ya tenemos el diplotipo de los 660 pacientes 
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP1A2_Allele_def_rs_excluidos)


### Pasos CYP3A4
Tras correr el bucle de script_genostar, se genera el objeto R CYP3A4_resultados_final. Hay info de diplotipo para 620 individuos y 40 no tiene diplotipo (NA)

Se sacan los genotipos de esos 40 individuos en df_CYP3a4_NA para poder trabajar con ellos 
#### Filtrar por los rs que usa genostar
Se filtra de los 40 individuos (df_CYP3A4_NA.csv) solo los rs que aparecern en la tabla de referencia que usa genostar, sacada de github genostar (CYP3A4_Allele_def.rda). Todo esta en el script filtrar_rs_referencia_genostar.
#### Correr el genostar de los NA: assign_diplotype 
Se ve que haciendo el filtrado se saca el diplotipo de los 40 pacientes restantes (CYP3A4_resultados_final_NA.rds), por loq eu ya tenemos el diplotipo de los 660 pacientes 
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP3A4_Allele_def_rs_excluidos)

### Pasos CYP3A5
Como hemos visto que filtrando para los rs que estan en la referencia de genostar si que funciona y va rapido, para este CYP se va a filtrar antes y luego correr genostar. Se hace todo en filtrar_rs_referencia_genostar
#### Filtrar por los rs que usa genostar
Se filtra de TODOS los individuos (df_CYP3A5.csv) solo los rs que aparecern en la tabla de referencia que usa genostar, sacada de github genostar (CYP3A5_Allele_def.rda). Todo esta en el script filtrar_rs_referencia_genostar.
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP3A5_Allele_def_rs_excluidos)
#### Correr el genostar de TODOS: assign_diplotype 
Se saca el diplotipo de loS 660 individuos (CYP3A5_resultados_final_diplotipos)

#### Correr el genostar de TODOS: all_geno_pheno 
Da un error al intentar sacar el metabolizador 

### Resultados
-Para CYP1A2 ya se tienen los diplotipos de 609 individuos en CYP1A2_resultados_final.rds y los 51 individuos restantes en CYP1A2_resultados_final_NA.rds, todos sacados con genostar 
-Para el CYP3A4 ya se tienen los diplotipos de  620 individuos en CYP3A4_resultados_final.rds y de los 40 individuos restantes en CYP3A4_resultados_final_NA.rds, todos sacados con genostar. 
-Para CYP3A5 ya se tienen los diplotipos de los 660 individuos en CYP3A5_resultados_final_diplotipos.rds, pero no los metabolizadores 



## 26/03/26
### Objetivo
Reunión con Claudio 
### Resultados
Me dijo que hasta el momento bien, pero que antes de seguir, me dijo de ir gen por gen y ver porque había rs que nosotros tenemos genotipados pero no usa genostar y ver porque genostar usa rs que no tenemos nosotros. Si es información redundante o si se esta dejando información importante fuera 


## 30/03/26
### Objetivo
Analizar rs de CYP1A2 y reunión con Claudio 
### Resumen rs CYP1A2
Hay 22 rs en comun entre genostar y los nuestros. Nosotros tenemos 60 más que genostar no usa. Genostar usa 24 más que nosotros no tenemos 
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP1A2_Allele_def_rs_excluidos)

### Resultados
Me dijo que hasta el momento bien, pero que antes de seguir, me dijo de ir gen por gen y ver porque había rs que nosotros tenemos genotipados pero no usa genostar y ver porque genostar usa rs que no tenemos nosotros. Si es información redundante o si se esta dejando información importante fuera 



## 31/03/26
### Objetivo
Hacer tabla CYP1A2_TABLA_COMPLETA  
### Resumen rs CYP1A2
-Hay 22 rs en comun entre genostar y los nuestros. Nosotros tenemos 60 más que genostar no usa. Genostar usa 24 más que nosotros no tenemos. En total hay 106 rs 
#### Tabla CYP1A2_TABLA_COMPLETA 
Se hace con BUSCARV 
-rsID: resultados de R 
-POS (hg19): se saca con BUSCARV a hoja SNP (se hace un liftover tambien)
-POS (hg38): se saca con BUSCARV a hoja SNV 

-Genostar: con BUSCARV de resulatdos de R 

-Nuestros datos: con BUSCARV de resultados de R

-1K_genomes: se saca con BUSCARV a hoja 1K_genomes (se saca de subset_g1000_chr15.txt al filtrar g1000_noFinnish.bim por cromosoma 15 y las posiciones del CYP1A2 de CYPSconventana19SNP)
```bash
awk '$1 ~ /^15$/ && $4 > 75036186 && $4 < 75053948 {print}' g1000_noFinnish.bim > subset_g1000_chr15.txt
```
-SNP: se saca con BUSCARV a hoja SNP (se saca filtrando por cromosoma 15 el 660indivCYPS_SNPs_19_02_26)

-SNV: se saca con BUSCARV a hoja SNV (se saca filtrando por cromosoma 15 el 660indivCYPS_SNVs_05_03_26)

-Genotipado: todas las SNVs estan imputadas y para las SNPs se sacan con BUSCARV a hoja genotipadas (se saca al filtrar Variantes_genotipas_SNPs por cromosoma 15). Undefined para las que no tenemos

-Imputado: al reves que genotipado. Undefined para las que no tenemos 

-MAF: se saca con BUSCARV a hoja MAF (del CYP1A2_rs_haploview_data_export, que se saca:)
##### 1. Crear CYP1A2_rs.txt con los 106 rs
##### 2. Sacar ped e info file con plink 
Se usan como bfile los g1000_noFinnish (sacado de marte) para solo las rs de CYP1A2_rs.txt 
```bash
plink --bfile g1000_noFinnish --extract CYP1A2_rs.txt --recode HV --out CYP1A2_haploview  
```
Solo se queda con las 61 rs que estan en 1K_genomes porque ha usado g1000_noFinnish
##### 3. Meter en haploview la info
En LINKAGE FORMAT se mete como data file: CYP1A2_haploview.chr-15.ped y como locus info: CYP1A2_haploview.chr-15.info
##### 4. Exportar la tabla de Check Markers como CYP1A2_rs_haploview_data_export

-HWpval: se saca con BUSCARV a hoja MAF 

-HOJA LD: se saca del CYP1A2_rs_LD.ld, que se saca:
##### 1. Crear CYP1A2_rs.txt con los 106 rs (ya esta hecho antes)
##### 2. Sacar ped e info file con plink 
Se usan como bfile los g1000_noFinnish (sacado de marte) para solo las rs de CYP1A2_rs.txt 
```bash
plink --bfile g1000_noFinnish --extract CYP1A2_rs.txt --r2 --out CYP1A2_rs_LD 
```
Solo se queda con las 61 rs que estan en 1K_genomes porque ha usado g1000_noFinnish
##### 3. Meter CYP1A2_rs_LD.ld en la tabla

### Resultados
Se termina la tabla CYP1A2_TABLA_COMPLETA


## 06/04/26
### Objetivo
Leer artículo PharmGKB CYP1A2 y analizar tabla CYP1A2_TABLA_COMPLETA  
### Artículo PharmGKB CYP1A2
"PharmGKB summary: very important pharmacogene information for CYP1A2" del 2012
Se ha visto que esta en nomenclatura old y hay tres rs importantes:
-rs762551: esta en nuestros datos yu en genostar 
-rs2069514: esta en nuestros datos pero no en genostar. Se intenta ver si se podría inferir con alguno de los otros rs de genostar: es el  rs9 y con haploview se ve que tiene r2 de 73 con los rs50 y rs 51 pero no estan en genostar, al igual que un r2 de 22 con el rs5 que tampoco esta en genostar, por lo que no sepuede inferir 
-rs12720461: esta en nuestros datos pero no en genostar. Se hizo lo mismo con el haploview y salian muchos r2 de 22 pero ninguno estaba en genostar 
### Resumen rs CYP1A2
-Hay 22 rs en comun entre genostar y los nuestros. Nosotros tenemos 60 más que genostar no usa. Genostar usa 24 más que nosotros no tenemos. En total hay 106 rs 
-La tabla CYP1A2_TABLA_COMPLETA ya esta hecha 
#### De los 24 rs que estan en genostar pero no en nuestros datos 
-Analizando la primera hoja hemos visto que solo hay tres variantes (marcadas en amarillo) que estan en genostar y en 1K_genomes pero no en nuestros datos. Al analizarlas en la hoja LD no salen por lo que no se puede inferir nada de ellas, son singleton 
-Analizando el excel CYP1A2_Allele_def_rs_excluidos (sacado con srcipt filtrar_rs_referencia_genostar), se hace una nueva hoja (filtrada) para ver que haplotipos estamos dejando fuera al no tener esos 24 rs en nuestros datos, son 24 haplotipos. Se sacan las frecuencias en la poblacion de estos 24 haplotipos en el excel de CYP1A2_Allele_def_rs_excluidos en la hoja de filtradas, como no hay info en ClinPGX se usa GenomAD y se usa el rs asociado a ese *alelo. Y se ve que todos los *alelos que no se tienen en cuenta tienen un frecuencia muy baja en la poblacion europea (non-Finnish). Hay uno que no se encuentra X3468AC porque tiene nomenclatura rara (marcado en rojo)
#### De los 60 rs que estan en nuestros datos pero no en genostar
-Hay 8 rs que no estan en 1K_genomes y no se puede hacer nada 
-rs2069514 y rs12720461 (del articulo del 2012), resaltados en rojo, estan en nuestros datos pero no en genostar, se intenta ver si se podrian inferir con alguno de los otros rs de genostar con el r2 del haploview, pero con lo que tienen algo de r2 tampoco estan en genostar, asi que son singleton 
-En la hoja LD todas las relaciones que salen, los rs estan en nuestros datos pero no en genostar, asi que no se puede hacer nada, son singletons  
### Resultados
Se terminan las limitaciones del CYP1A2 



## 07/04/26
### Objetivo
Hacer y analizar tabla CYP3A4_TABLA_COMPLETA  
### Resumen rs CYP3A4
-Hay 18 rs en comun entre genostar y los nuestros. Nosotros tenemos 70 más que genostar no usa. Genostar usa 24 más que nosotros no tenemos. En total hay 112 rs 
-Se hace la tabla CYP3A4_TABLA_COMPLETA como para CYP1A2
#### De los 24 rs que estan en genostar pero no en nuestros datos 
-Analizando la primera hoja hemos visto que solo hay dos variantes (marcadas en amarillo) que estan en genostar y en 1K_genomes pero no en nuestros datos. Al analizarlas en la hoja LD no salen por lo que no se puede inferir nada de ellas, son singleton 
-Analizando el excel CYP3A4_Allele_def_rs_excluidos (sacado con srcipt filtrar_rs_referencia_genostar), se hace una nueva hoja (filtrada) para ver que haplotipos estamos dejando fuera al no tener esos 24 rs en nuestros datos, son 24 haplotipos. Se sacan las frecuencias en la poblacion de estos 24 haplotipos en el excel de CYP3A4_Allele_def_rs_excluidos en la hoja de filtradas, como no hay info en ClinPGX se usa GenomAD y se usa el rs asociado a ese *alelo. Y se ve que todos los *alelos que no se tienen en cuenta tienen un frecuencia muy baja en la poblacion europea (non-Finnish), hay algunos que se buscan en el UCSC (marcados en amarillo) y cinco que no se encuentran, cuatro porque tienen nomenclatura rara (marcados en rojo) y uno porque no esta su frecuencia en ninguna base de datos (marcado en morado)
#### De los 60 rs que estan en nuestros datos pero no en genostar
-Hay 15 rs que no estan en 1K_genomes y no se puede hacer nada 
-En la hoja LD todas las relaciones que salen, los rs estan en nuestros datos pero no en genostar, asi que no se puede hacer nada, son singletons  
### Resultados
Se terminan las limitaciones del CYP3A4 



## 08/04/26
### Objetivo
Hacer repositorio Github en el ordenador de aqui y en el mio 
### Resultados
Se hace el repositorio bien 



## 09/04/26
### Objetivo
Mirar las limitaciones de CYP3A5 y sacar diplotipos y limitaciones de CYP2C19, CYP2C9, CYP2B6 y CYP2D6
### Resumen rs CYP3A5
-Hay 4 rs en comun entre genostar y los nuestros. Nosotros tenemos 73 más que genostar no usa. Genostar usa 1 más que nosotros no tenemos. En total hay 78 rs 
-Ya no se hace la TABLA COMPLETA para CYP3A5 porque se ha visto que no da info 
#### Del rs que esta en genostar pero no en nuestros datos (limitaciones)
-Analizando el excel CYP3A5_Allele_def_rs_excluidos (sacado con srcipt filtrar_rs_referencia_genostar), se hace una nueva hoja (filtrada) para ver que haplotipos estamos dejando fuera al no tener ese rs en nuestros datos, es 1 haplotipo. Se sacan la frecuencia en la poblacion de ese haplotipos en el excel de CYP3A5_Allele_def_rs_excluidos en la hoja de filtradas, se saca de CYP3A5_frequency_table (descaragdo de ClinPGX pharmacogene table frequencies), se mira ese *alelo para Europeos y se ve que no tiene frecuencia asi que bien. 


### Pasos CYP2B6
En filtrar_rs_referencia_genostar: primero se filtra por rs comunes para que sea más rápido (se guarda CYP2B6_Allele_def_rs_excluidos para las limitaciones) y luego se aplica genostar assign_diplotype (se guarda en CYP2B6_resultados_final_diplotipos). Al usar all_geno_pheno para metabolizador da error
### Resumen rs CYP2B6
-Hay 30 rs en comun entre genostar y los nuestros. Nosotros tenemos 232 más que genostar no usa. Genostar usa 18 más que nosotros no tenemos. En total hay 280 rs 
-Ya no se hace la TABLA COMPLETA para CYP2B6 porque se ha visto que no da info 
#### De los 18 que esta en genostar pero no en nuestros datos (limitaciones)
-Analizando el excel CYP2B6_Allele_def_rs_excluidos (sacado con srcipt filtrar_rs_referencia_genostar), se hace una nueva hoja (filtrada) para ver que haplotipos estamos dejando fuera al no tener esos rs en nuestros datos, son 21 haplotipos. Se sacan la frecuencia en la poblacion de esos 21 haplotipos en el excel de CYP2B6_Allele_def_rs_excluidos en la hoja de filtradas, se saca de CYP2B6_frequency_table cinco frecuencias (descaragdo de ClinPGX pharmacogene table frequencies) las marcadas en verde, se miran esos *alelos para Europeos. Luego hay otras que no tienen info en ClinPGX se usa GenomAD y se usa el rs asociado a ese *alelo, hay cuatro que no se encuentran porque tienen nomenclatura rara (marcados en rojo). El error es muy alto 
#### Diplotipos CYP2B6
En CYP2B6_resultados_final_diplotipos


### Pasos CYP2C9
En filtrar_rs_referencia_genostar: primero se filtra por rs comunes para que sea más rápido (se guarda CYP2C9_Allele_def_rs_excluidos para las limitaciones) y luego se aplica genostar assign_diplotype (se guarda en CYP2C9_resultados_final_diplotipos). Al usar all_geno_pheno para metabolizador da error
### Resumen rs CYP2C9
-Hay 41 rs en comun entre genostar y los nuestros. Nosotros tenemos 262 más que genostar no usa. Genostar usa 39 más que nosotros no tenemos. En total hay 342 rs 
-Ya no se hace la TABLA COMPLETA para CYP2C9 porque se ha visto que no da info 
#### De los 39 que esta en genostar pero no en nuestros datos (limitaciones)
-Analizando el excel CYP2C9_Allele_def_rs_excluidos (sacado con srcipt filtrar_rs_referencia_genostar), se hace una nueva hoja (filtrada) para ver que haplotipos estamos dejando fuera al no tener esos rs en nuestros datos, son 40 haplotipos. Se sacan la frecuencia en la poblacion de esos 40 haplotipos en el excel de CYP2C9_Allele_def_rs_excluidos en la hoja de filtradas, se saca de CYP2C9_frequency_table ocho frecuencias (descaragdo de ClinPGX pharmacogene table frequencies) las marcadas en verde, se miran esos *alelos para Europeos. Luego hay otras que no tienen info en ClinPGX se usa GenomAD y se usa el rs asociado a ese *alelo, hay 18 que no se encuentran porque tienen nomenclatura rara (marcados en rojo) 
#### Diplotipos CYP2C9
En CYP2C9_resultados_final_diplotipos



### Pasos CYP2C19
En filtrar_rs_referencia_genostar: primero se filtra por rs comunes para que sea más rápido (se guarda CYP2C19_Allele_def_rs_excluidos para las limitaciones) y luego se aplica genostar assign_diplotype (se guarda en CYP2C19_resultados_final_diplotipos). Al usar all_geno_pheno para metabolizador da error
### Resumen rs CYP2C19
-Hay 24 rs en comun entre genostar y los nuestros. Nosotros tenemos 303 más que genostar no usa. Genostar usa 11 más que nosotros no tenemos. En total hay 338 rs 
-Ya no se hace la TABLA COMPLETA para CYP2C19 porque se ha visto que no da info 
#### De los 11 que esta en genostar pero no en nuestros datos (limitaciones)
-Analizando el excel CYP2C19_Allele_def_rs_excluidos (sacado con srcipt filtrar_rs_referencia_genostar), se hace una nueva hoja (filtrada) para ver que haplotipos estamos dejando fuera al no tener esos rs en nuestros datos, son 12 haplotipos. Se sacan la frecuencia en la poblacion de esos 12 haplotipos en el excel de CYP2C19_Allele_def_rs_excluidos en la hoja de filtradas, se saca de CYP2C19_frequency_table dos frecuencias (descaragdo de ClinPGX pharmacogene table frequencies) las marcadas en verde, se miran esos *alelos para Europeos. Luego hay otras que no tienen info en ClinPGX se usa GenomAD y se usa el rs asociado a ese *alelo, hay cinco que no se encuentran, dos porque tienen nomenclatura rara (marcados en rojo) y tres porque no esta su frecuencia en ninguna base de datos (marcado en morado) 
#### Diplotipos CYP2C19
No salen todo, salen muchos unknown (CYP2C19_resultados_final_diplotipos), deberíamos tener en cuenta los CNVs 



### Pasos CYP2D6
En filtrar_rs_referencia_genostar: primero se filtra por rs comunes para que sea más rápido (se guarda CYP2D6_Allele_def_rs_excluidos para las limitaciones) y luego se aplica genostar assign_diplotype (se guarda en CYP2D6_resultados_final_diplotipos) pero no sale creo que es por las CNVs 
### Resumen rs CYP2D6
-Hay 40 rs en comun entre genostar y los nuestros. Nosotros tenemos 77 más que genostar no usa. Genostar usa 110 más que nosotros no tenemos. En total hay 227 rs 
-Ya no se hace la TABLA COMPLETA para CYP2B6 porque se ha visto que no da info 
-Necesita CNVs
-No se calcula el error, hasta meter las CNVs 


### Resultados
-Se terminan las limitaciones del CYP3A5 y CYP2B6 bien. El error de CYP2B6 es muy alto mirar la razon 
-Tambien se hace el error para CYP2C19 y CYP2C9
-Para CYP2C19 los diplotipos salen mal, deberían mirarse las CNVs
-Para CYP2D6 no sale nada, se necesitan los CNVs 



## 10/04/26
### Objetivo 
-CNVs 
-Organizar carpetas 
### CNVs en CNV_gene_overlap_database en la hoja filtered (SON PERSONAS)
Hat 125 en total en filtradas 
-CYP1A2: nada 
-CYP3A4: 1 duplicacon 
-CYP3A5: 1 duplicacion 
-CYP2C19: 19 duplicaciones y 2 deleciones
-CYP2C9: 19 duplicaciones
-CYP2B6: 1 duplicacion y 1 delecion 
-CYP2D6: 79 duplicaciones y 2 deleciones 
### Resultados
Hay muchos individuos con CNVs en CYP2D6, CYP2C19 y CYP2C9, hay que mirar si afecta en genostar o no 



## 13/04/26
### Objetivo  
-Tabla general resumen 
-Mirar el error alto de CYP2B6  
### Tabla general resumen 
Se juntan todos los datos de los CYPS en TABLA_RESUMEN_CYPS 
-rs total, rs comunes, rs solo nosotros, rs solo genostar: se saca de la info apuntada en este documento
-nº *alelos no tenidos en cuenta en genostar: se saca de CYPX_Allele_def_rs_excluidos. Se añade esta columna porque un rs (de solo genostar) puede dar a más de un alelo
-nº *alelos sin frec poblacional de los *alelos no tenidos en cuenta en genostar: se saca de CYPX_Allele_def_rs_excluidos (marcados en rojo o morado)
-Error rs solo genostar : se saca de CYPX_Allele_def_rs_excluidos
### Mirar el error alto de CYP2B6
Un 14.02% del 14.07% del error se debe a la variante rs3211371 que define el *5 y *7 muy frecuentes en la población europea, se cree que es que no esta genotipada. Al buscar en Variantes_genotipas_SNPs no sale y en concatenar_CYPS, que son todas las variantes que tenemos en cuenta. Se me ocurre mirar si esta en LD alto con alguna de las que si que tenemos genotipadas. 
#### MAF y LD y haploview 
Se hace la tabla CYP2B6_TABLA_ERROR_LD
-rsID: resultados de R 
-Genostar: con BUSCARV de resulatdos de R 

-Nuestros datos: con BUSCARV de resultados de R

-1K_genomes: se saca con BUSCARV a hoja 1K_genomes (se saca de subset_g1000_chr19.txt al filtrar g1000_noFinnish.bim por cromosoma 19 y las posiciones del CYP2B6 de CYPSconventana19SNP)
```bash
awk '$1 ~ /^19$/ && $4 > 41492187 && $4 < 41529303 {print}' g1000_noFinnish.bim > subset_g1000_chr19.txt
```
Nos salen 211 rs que estan tanto en nuestros datos como en 1K_genomes 

-MAF: se saca con BUSCARV a hoja MAF (del CYP2B6_rs_haploview_data_export, que se saca:)
##### 1. Crear CYP2B6_rs.txt con los 262 rs que nosotros tenemos y el rs3211371 
##### 2. Sacar ped e info file con plink 
Se usan como bfile los g1000_noFinnish (sacado de marte) para solo las rs de CYP2B6_rs.txt 
```bash
plink --bfile g1000_noFinnish --extract CYP2B6_rs.txt --recode HV --out CYP2B6_haploview  
```
Solo se queda con 212 rs, las 211 rs que estan en 1K_genomes porque ha usado g1000_noFinnish y la rs3211371
##### 3. Meter en haploview la info
En LINKAGE FORMAT se mete como data file: CYP2B6_haploview.chr-15.ped y como locus info: CYP2B6_haploview.chr-15.info
##### 4. Exportar la tabla de Check Markers como CYP2B6_rs_haploview_data_export

-HWpval: se saca con BUSCARV a hoja MAF 

-HOJA LD: se saca del CYP2B6_rs_LD.ld, que se saca:
##### 1. Crear CYP2B6_rs.txt con los 262 rs que nosotros tenemos y el rs3211371 
##### 2. Sacar ped e info file con plink 
Se usan como bfile los g1000_noFinnish (sacado de marte) para solo las rs de CYP2B6_rs.txt 
```bash
plink --bfile g1000_noFinnish --extract CYP2B6_rs.txt --r2 --out CYP2B6_rs_LD 
```
Solo se queda con 212 rs, las 211 rs que estan en 1K_genomes porque ha usado g1000_noFinnish y la rs3211371
##### 3. Meter CYP2B6_rs_LD.ld en la tabla

Al mirar el haploview es la 170 pero no tiene LD alto con ninguna de las nuestras y en la tabla de LD no aparece

#### Frecuencia población europea *5 y rs3211371 
Es la misma frec asi que perfecto 
-5*: 0.11547839 (CYP2B6_frequency_table)
-rs3211371: en GnomAD 0.1230 (https://gnomad.broadinstitute.org/variant/19-41016810-C-T?dataset=gnomad_r4), en 1K_genomes (hoja MAF de CYP2B6_TABLA_ERROR_LD) 0.1230, en dbSNP 0.09212 (https://www.ncbi.nlm.nih.gov/snp/rs3211371)

### Resultados 
Se hace la TABLA_RESUMEN_CYPS (falta CYP2D6 pero cuando tengamos las CNVs nuevas) y todos los errores bien excepto para CYP2B6, que la culpable es rs3211371. Info sobre esta variante:
-La usa genostar pero nosotros no la tenemos 
-Si que se genotipó pero por filtrado se eliminó
-Es Trialelica 
-Da el *5 sola, *7 con otras variantes y *33 y *34 (pero no contribuyen al error)
-El *5 es muy frecuente en poblacion europea pero da METABOLIZADOR NORMAL
-La frecuencia de *5 y de rs3211371 es la misma 
Se llega a la conclusion de que aunque la frecuencia es alta, da metabolizador normal por lo que no importa mucho 


## 14/04/26
### Objetivo
-Sacar frecuencias diplotipos, allele star para CYP1A2, CYP3A4, CYP3A5, CYP2C19, CYP2C9 y CYP2B6 y metabolizadores para CYP3A5, CYP2C19, CYP2C9 y CYP2B6
### Se generan tablas y gráficas 
#### Tabla CYP_tabla_diplotipos: ID Y DIPLOTIPO 
Con cualquiera de los scripts, usa de input CYP_resultados_final_diplotipos 
#### Plot CYP_plot_frec_diplotipos: FRECUENCIA ABSOLUTA DIPLOTIPOS
Con script_plot_frec_diplotipos, usa de input CYP_resultados_final_diplotipos
Se ven todos los diplotipos grises
#### Plot CYP_plot_prop_diplotipos: PROPORCIONES DIPLOTIPOS  
Con script_plot_prop_diplotipos, usa de input CYP_resultados_final_diplotipos
Se ven en rojo los diplotipos NA, en verde el *1/*1 y en azul otros diplotipos 
#### Plot CYP_plot_prop_star_allele: PROPORCIONES STAR ALLELES  
Con script_plot_prop_allele_star, usa de input CYP_resultados_final_diplotipos
Se ven en rojo los NA, en verde el *1 y en azul otros *alelos
#### Plot CYP_plot_prop_diplotipos_metabolizador: PROPORCIONES DIPLOTIPOS coloreados por METABOLIZADOR 
Con script_plot_prop_diplotipos_metabolizadores, usa de input CYP_resultados_final_diplotipos y las CYP_Diplotype_Phenotype_Table (pharmacogene tables de ClinPGX)
Se colorean los diplotipos según el metabolizador  
#### Plot CYPs_plot_prop_metabolizador: PROPORCIONES METABOLIZADOR todos los CYPs 
Con script_plot_prop_metabolizador_CYPS, usa de input CYP_resultados_final_diplotipos y las CYP_Diplotype_Phenotype_Table (pharmacogene tables de ClinPGX)
### Resultados
Se saca 
-CYP_tabla_diplotipos: ID Y DIPLOTIPO
-CYP_plot_frec_diplotipos: FRECUENCIA ABSOLUTA DIPLOTIPOS
-CYP_plot_prop_diplotipos: PROPORCIONES DIPLOTIPOS
-CYP_plot_prop_star_allele: PROPORCIONES STAR ALLELES 
-CYP_plot_prop_diplotipos_metabolizador: PROPORCIONES DIPLOTIPOS coloreados por METABOLIZADOR (no se hace para CYP1A2 y CYP3A4 porque no hay CYP_Diplotype_Phenotype_Table)
-CYPs_plot_prop_metabolizador: PROPORCIONES METABOLIZADOR todos los CYPs 

Estaria bien mirar los metabolizadores de CYP1A2 y CYP3A4


## 15/04/26
### Objetivo
-Charla prevención 
-Comparar si las proporciones que nos han salido de *allele star se corresponden con las frecuencias poblacionales esperadas 
### CYP3A5
Se compara CYP3A5_plot_prop_star_allele con CYP3A5_frequency_table y cuadra perfecto 
### CYP2C9
Se comparan CYP2C9_plot_prop_star_allele con CYP2C9_frequency_table y cuadra perfecto 
### CYP2B6
Se comparan CYP2B6_plot_prop_star_allele con CYP2B6_frequency_table, cuadra más o menos excepto el *5 que no esta en nuestros datos y es frecuente en la población 
### CYP1A2
No hay tabla referencia. CYP1A2_Allele_def_2.rda se pasa a CYP1A2_Allele_def_rs_star_allele, que se filtra por los *alelos que hay en CYP1A2_plot_prop_star_allele y se mira la frecuencia de los rs que los definen. Y cuadran bastante bien con las frecuencias poblacionales 
### CYP3A4
No hay tabla referencia. CYP3A4_Allele_def.rda se pasa a CYP3A4_Allele_def_rs_star_allele, que se filtra por los *alelos que hay en CYP3A4_plot_prop_star_allele y se mira la frecuencia de los rs que los definen. Y cuadran bastante bien con las frecuencias poblacionales 
### CYP2C19
No se sacaron bien los diplotipos de este CYP, mirar con CNVs 
### Resultados
Las frecuencias de los allele star que nos han salido se corresponden con lo esperado en la poblacion bastante bien  


## 16/04/26
### Objetivo 
-CNVs nuevo 
-Reunion Claudio 
-Mirar los rs raros rojos ver si los tenemos por posicion 
### CNVs en CNV_gene_overlap_database_nuevo (SON PERSONAS)
Hat 637 en total en filtradas 
-CYP1A2: 2 duplicaciones 
-CYP3A4: 99 duplicaciones 
-CYP3A5: 78 duplicaciones 
-CYP2C19: 89 duplicaciones y 2 deleciones
-CYP2C9: 78 duplicaciones
-CYP2B6: 252 duplicacion y 1 delecion 
-CYP2D6: 34 duplicaciones y 2 deleciones 
### Ver por posicion si tenemos los rs raros 
1. Comparar CYP_Allele_def.rda (la que usa genostar) con CYP_allele_definition_table (ClinPGX) y ver la posicion de los X (rojos) de CYP_Allele_def_rs_excluidos
2. Buscar por posición si nosotros tenemos ese rs en 660indivCYPS_SNPs_SNVs_05_03_26.bim
#### CYP1A2 
Hay 1 rojo en CYP_Allele_def_rs_excluidos, Se compara CYP1A2_Allele_def_2.rda con CYP1A2_Haplotypes-PS216394-1454147960 y se ve que la que estaba en rojo correponde al rs2505559894, que no esta en nuestros datos, por lo que se mira su frecuencia y se pone en azul. Frecuencias completas   
#### CYP3A5
Se compara CYP3A5_Allele_def.rda con CYP3A5_allele_definition_table y se ve que todo bien, que no hay cosas raras. Frecuencias completas 
#### CYP2C9
Hay varias rojo en CYP2C9_Allele_def_rs_excluidos, Se compara CYP2C9_Allele_def.rda con CYP2C9_allele_definition_table y se ve que varias en rojo correponden a rs que no estan en nuestros datos, por lo que se mira su frecuencia y se pone en azul, luego hay otra que se saca por posicion pero no la tenemos y tambien se pone en azul, luego hay otras que no tienen rs asociado y se mira por pisicion en nuestros datos y tampoco asi que se quedan en rojo, pero estan en ClinPGX asi que se pone la frecuencia del *alelo. Frecuencias completas 
#### CYP2C19
Hay varias rojo en CYP2C19_Allele_def_rs_excluidos, Se compara CYP2C19_Allele_def.rda con CYP2C19_allele_definition_table y se ve que algunas en morado y rojo estan en ClinPGX asi que se pone la frecuencia del *alelo. Frecuencias completas 
#### CYP2B6
Hay varias rojo en CYP2B6_Allele_def_rs_excluidos, Se compara CYP2B6_Allele_def.rda con CYP2B6_allele_definition_table y se ve que una en rojo correponde a rs que no estan en nuestros datos, estan en ClinPGX asi que se pone la frecuencia del *alelo. Luego hay otra (X) que no tiene rs asociado y se mira por pisicion en nuestros datos y tampoco asi que se quedan en rojo. Nos queda un rs sin frecuencia (X)
#### CYP3A4
Hay 4 rojo en CYP_Allele_def_rs_excluidos, Se compara CYP3A4_Allele_def.rda con CYP3A4_allele_definition_table-PS216409-1454596980 y se ve que dos en rojo correponden a rs que no estan en nuestros datos, por lo que se mira su frecuencia y se pone en azul, luego hay dos que no tienen rs asociado y se mira por pisicion en nuestros datos y tampoco asi que se quedan en rojo. Algunas estan en ClinPGX asi que se pone la frecuencia del *alelo. Nos quedan 3 rs sin frecuencia (X.3, X.2 y rs1318364992)
### Resultados
Hay muchos individuos con CNVs en CYP2D6, CYP2C19 y CYP2C9, hay que mirar si afecta en genostar o no 



## 17/04/26
### Objetivo 
Organizar un poco todo 
### Resumen CYPs del error 
-CYPA12: todas las frecuencias de los rs excluidos sacadas   
-CYP2C9: todas las frecuencias de los rs excluidos sacadas  
-CYP3A5: todas las frecuencias de los rs excluidos sacadas 
-CYP3A4: hay tres rs de los excluidos que no se puede sacar su frecuencia 
-CYP2B6: hay un rs de los excluidos que  no se puede sacar su frecuencia, tiene VARIANTES ESTRUCTURALES (X.2)
-CYP2C19: todas las frecuencias de los rs excluidos sacadas, tiene VARIANTES ESTRUCTURALES (X) y cuidado que el * de ref es el 38 pero genostar lo cambia 
-CYP2D6: mirar los CNVs 
### Resulatdos 
Mirar las variantes estructurales 

## 20/04/26
### Objetivo 
-Jugar con genostar para que me saque lo máximo posible en una sola función 
#### script_juntar_dfs_CYPS: genera matriz_final_sin_CNVs_CYP2D6
Coge todas las tablas de referencia que usa genostar para cada CYP y matrix_geno_fixed_espacios y saca los rs comunes en nuestros datos y genostar para cada CYP y los junta todos con el formato correcto para genostar 
#### script_prueba_funciones_genostar y script_todo_genostar: saca los diplotipos, haplotipos, metabolizadores y pie charts 



## 21/04/26
### Objetivo 
-Apuntar dudas genostar 
-Leer información Genostar del github (las funciones de man y en R el genotype_conversion y en data el allele_definitions_snapshot_2025-08-13) 
### Resultados 
-Se hace una jerarquía de las funciones y se hace un resumen general de cada una (en el word Funciones_genostar)
-Apuntar dudas en el cuaderno 


## 22/04/26
### Objetivo 
-Gestionar las CNVs de CYP2D6
-Probar genostar con la matriz 
### CNVs de CYP2D6
Se hacen a mano con el custom track de UCSC (CYP2D6_CNVs), se hace para CNVX9 (exon 9), CNVInt6 (intron 6) y CNV5UTR (zona 5UTR) y solo nos salen 3 individuos con CNVs de los 660 totales: 
-HSP_0041.CEL: duplicacion en CNVX9 (exon 9), CNVInt6 (intron 6) y CNV5UTR (zona 5UTR)
-PNT_0196.CEL: duplicacion en CNVX9 (exon 9), CNVInt6 (intron 6) y CNV5UTR (zona 5UTR)
-UAM_0343.CEL: delecion en CNVX9 (exon 9) y CNVInt6 (intron 6), normal en CNV5UTR (zona 5UTR)
### Juntar CYP2D6_CNVs con matriz_final_sin_CNVs_CYP2D6
En el script script_juntar_dfs_CYPS 
### Resultados 
-Se deja la matriz MATRIZ_FINAL_COMPLETA ya hecha con todos los CYPs y la CNVs de CYP2D6 incluida
-Se prueba genostar con la matriz pero salen errores
1. Al correr CYP2D6 entero 
lista_diplotipos <- all_geno_pheno(df, c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6"), phased = FALSE, CYP1A2_name = "new")
Error in rs1065852 == "AA" | "A A" : solo son posibles operaciones para variables de tipo numérico, compleja o lógico
SE PRUEBA A QUITAR ESE RS  
2. Al correr all_geno_pheno 
lista_diplotipos <- all_geno_pheno(df_sin_rs, c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6"), phased = FALSE, CYP1A2_name = "new")
Error del split
SE HACE ASSIGN_DIPLOTYPE
3. Correr assign_diplotype
lista_diplotipos <- assign_diplotype(df_sin_rs, c("CYP3A5","CYP1A2","CYP3A4","CYP2C9","CYP2B6","CYP2C19","CYP2D6"), phased = FALSE, CYP1A2_name = "new")
SALE BIEN, PERO PARA CYP2C19 Y CYP2D6 MUCHOS NA 


## 23/04/26
### Objetivo
-Reunión Claudio 
-Sacar cosas con genostar
-Corregir los errores de los CYPS   
-Mirar los problemas del día anterior 
### Sacar cosas con genostar 
En script_genostar  
### Corregir los errores de los CYPS 
En CYP_Allele_def_rs_excluidos, las frecuencias que se han sacado por *alelo hay que verificar que no haya mas rs que dan a ese *alelo o que la frecuencia del rs coindica con la frecuencia del *alelo. Se cambian estas tablas y los errores 
### Mirar los problemas del día anterior
Se ponen las funciones de genostar en funciones_genostar_diana y se cambian 
-adjust_diplotype: se cambia rs1065852 == "AA" | "A A" porque esta mal 
### Resultados 
Se actualizan los errores 




## 27/04/26
### Objetivo
-Actualizar errores TABLA_RESUMEN y actualizar error CYP2B6 
-Error de CYP2D6 
### Sacar cosas con genostar 
En script_genostar se saca todo con all_geno_pheno y se gurada todo el resultado en salida_genostar y el resultado filtrado con solo las columnas que interesan en salida_genostar_filtrado
### Error CYP2D6  
En CYP2D6_Allele_def_rs_excluidos. Hay rs problematicos anotados en el cuaderno   
### Resultados 
Se hace el error de CYP2D6 y se mete en TABLA_RESUMEN_CYPS


## 11/05/26
### Objetivo 
Hacer script de la grafica final de frecuencias (ClinPGX vs nuestros resultados)  
### Datos poblacionales: Excel frecuencias_poblacionales_metabolizadores_clinpgx.xlsx
Se hace un excel para CYP2C9, CYP2C19, CYP3A5, CYP2B6 y CYP2D6 con las frecuencias poblacionales de clinpgx (se saca de la pharmacogene table CYP_frequency_table.xlsx de la hoja Phenotype Frequency, para las que tienen activity score (CYP2C9 y CYP2D6) se usa CYP_Diplotype_Phenotype_Table para pasar de activity score a metabolizador)
### Nuestros datos: Excel salida_genostar_filtrado
Es la salida de script_genostar 
### Script_grafica_final_metabolizadores  
Se le mete 
-Las frecuencias de clinpgx (frecuencias_poblacionales_metabolizadores_clinpgx.xlsx) 
-Nuestra salida de genostar (salida_genostar_filtrado) 
Y te saca el plot resumen de frecuencias de metabolizadores (plot_final_metabolizadores) 
### Resultados
La mayoria cuadra, hay que revisar que los porcentajes esten bien y ver que pasa con CYP2C19 y CYP2D6 


### 12/05/26
### Objetivo
-Revisar plot_final_metabolizadores
-Hacer plot_final_allele_star
### Script_grafica_final_star_alleles  
Se le mete:
-Las frecuencias de clinpgx para CYP3A5, CYP2B6, CYP2C19 CYP2C9 y CYP2D6 (se saca de la pharmacogene table CYP_frequency_table.xlsx de la hoja Allele Frequency, se hace la columna European100 para expresar la frecuencia en porcentaje), mietras que para CYP1A2 y CYP3A4 no hay frecuencias en clinpgx se usan las frecuencias de los rs de esos allele star que se sacan de CYP1A2_Allele_def_rs_star_allele y CYP3A4_Allele_def_rs_star_allele y se ponen en frecuencias_poblacionales_allele_star_CYP1A2_CYP3A4_clinpgx
-Nuestra salida de genostar (salida_genostar_filtrado) 
Y te saca el plot resumen de frecuencias de allele star (plot_final_allele_star)  
### Resultados 
-El plot_final_metabolizadores esta bien, solo falta ver porque para CYP2C19 y CYP2D6 hay muchos NA 
-Se hace el plot_final_allele_star y esta bien 


### 13/05/26
### Objetivo
-Hacer plot_final_activity 
-Actualizar tabla TABLA_RESUMEN_CYPS
### Script_grafica_final_activity  
Se le mete 
-Las frecuencias de clinpgx para CYP2C9 y CYP2D6 (frecuencias_poblacionales_activity_clinpgx.xlsx: excel para CYP2C9 y CYP2D6 con las frecuencias poblacionales de clinpgx (se saca de la pharmacogene table CYP_frequency_table.xlsx de la hoja Phenotype Frequency))  
-Nuestra salida de genostar (salida_genostar_filtrado) 
Y te saca el plot resumen de frecuencias de activity (plot_final_activity) 
### Actualizar tabla TABLA_RESUMEN_CYPS con el CALL RATE 
Se calcula el call rate (porcentaje de individuos para los que no se consigue sacar diplotipo) de cada CYP con el script_call_rate 
### Resultados 
-Se hace el plot_final_activity_score y esta bien 
-Se actualiza TABLA_RESUMEN_CYPS con el call rate de cada CYP 




### 14/05/26
### Objetivo
-Analizar funciones en busca de la razon del bajo call rate para CYP2C19 y CYP2D6 
### Mejora en la funcion filter_no_match 
Se le añade para aparte de "No matching star alleles found" tambien funcione si hay NA (reales o "NA") y se sacan los excels con los individuos a los que no se les ha encontrado diplotipo CYP_individuos_NA
### Modificar las tablas de individuos NA tras genostar CYP_individuos_NA
-En la hoja salida genostar: se deja tal cual sale (con todas las columnas)
-En la hoja filtrada: se filtra para los rs y columnas de ese CYP 
### Buscar la razon por la que salen los NA 
Se deja en el cuaderno marcado pero en resumen en CYP_individuos_NA 
-Blanco: NA a los que le falta al menos un genotipo 
-Naranja: NA pero que tienen todos los genotipos 
### Resultados 
-Se hacen mejoras de las fucniones (en el cuaderno en la parte de mejoras de las funciones genostar)
-Para CYP3A4 y CYP2C9 los NA que salen se deben a que les falta algun dato de genotipado par al menos un rs 
-Para CYP2B6 la mayoria de los NA es porque no tienen algun dato de genotipado, excepto 2 que no se sabe de donde viene porque tienen todos los genotipos 
-Para CYP2D6 y CYP2C19 pocos son NA por falta de algun genotipo pero la mayoria (433 para CYP2C19 y 311 para CYP2D6) no se sabe de donde viene porque tienen todos los genotipos 


### 18/05/26
### Objetivo
-Hacer abstracts para CBM y Reus
-Posiciones epi
-Ver que individuos de epi coinciden con los 660 individuos 
### Posiciones de los CYPs para epi 
CYP, cromosoma, start y end de UCSC Genome Browser (a start y end se le pone una VENTANA de 100.000) en GRCh38/hg38, id del gen de https://ensembl.org/ 
Resultado: CYPSconventana38epi se ponen las GRCh38/hg38 (se saca de la hoja hg38 ventana 100.000 (epi) de CYPStodo)
### Individuos comunes Genostar y epi (individuos_EPIC)
En genostar hay 660 individuos y en epi 238. Todos los de epi tienen genostar
### Calcular call rate de genostar para los 238 individuos de epi
#### 1. Filtrar salida genostar para los 238 individuos de epi (salida_genostar_filtrado_epi)
Se coge salida_genostar_filtrado y se filtran los 238 individuos con info de epi 
#### 2. Calular call rate de esos individuos (genostar_CALL_RATE_individuos_epi)
Se usa el script script_CALL_RATE sobre el excel salida_genostar_filtrado_epi
### Resultados 
-Se suben los abstracts a marte 
-Se le pasa a Ines las posiciones de los CYPs con ventana de 100000 (CYPSconventana38epi)
-Los 238 individuos de epi se han corrido en genostar
-Call rate de genostar para los 238 individuos de epi (genostar_CALL_RATE_individuos_epi)

### 19/05/26
### Objetivo
-Analizar funciones en busca de la razon del bajo call rate para CYP2C19 y CYP2D6  
### Mejoras en la funcion find_missing_data_pre 
-Se le añade para aparte de " " tambien funcione si hay NA (reales o "NA") o celda vacia "" 
-Se modofica para que acepte cualquier nombre como primera columna y no solo patient_index 
Resultado: datos_faltantes
CUIDADO: se hace otra fucnion find_missing_data_pre para aplicarla antes y poder sacar los valores faltantes, pero no se aplica en el pipeline de genostar porque da errores 
### Solución del bajo call rate de CYP2C19 
La MATRIZ_FINAL_COMPLETA se pasa a MATRIZ_FINAL_COMPLETA_CYP2C19 que es la que se usa, se modifica el 
-CYP2C19_rs17882687: de CC a AA 
-CYP2C19_rs113934938: de AA a GG
### Se vuelven a correr los scripts para actualizar 
### Resultados 
-Se hacen mejoras de las fucniones (en el cuaderno en la parte de mejoras de las funciones genostar)



### 20/05/26
### Objetivo
-Transcriptomica: sacar media y desviación por CYP y ver donde se expresa cada CYP 
-Epi: anotar las CpG que han salido segun promotor o isla CpG 
-Actualizar TFM 
### Transcriptómica raw: CYP_gene_counts_238_individuos (son gene_counts)
### Epi raw: CpGs_por_CYP_ventana100kb_con_betas (beta values)
### Actualizar TFM 
En la primer apágina del cuaderno poner lo que se hace para transcriptomica y epi 
### Transcriptomica: media y desviacion  
Los datos de transcriptómica raw son: 
Con el script_CYP_media_ds se calcula la media y la desviacion tipica de cada CYP y se guarda en CYP_gene_counts_238_individuos en la hoja CYP_media_ds y sale que CYP2D6 y CYP3A5 tienen algunos counts
### Transcriptomica: expresion de los CYPs 
Se busca en The Human Protein Atlas en Tissue:RNA los CYPs y se ve que todos se expresan en el hígado y ninguno en sangre 
-CYP2C19, CYP3A4, CYP2C9, CYP1A2 y CYP2B6: nada de expresión en blood 
-CYP3A5: tiene un poco de expresión en blood vessel https://www.proteinatlas.org/ENSG00000106258-CYP3A5/tissue 
-CYP2D6: NO tiene espresión en blood https://www.proteinatlas.org/ENSG00000100197-CYP2D6/tissue
### Epi: anotación de las CpG  
Se guarda en CpGs_por_CYP_ventana100kb_con_betas en sus ventanas correspondientes (se hacen las hojas cosas_CYP para resumir lo que se ve en la ventana)
En la hoja resumen CpGs se apunta las CpGs por gen 
### Resultados 
-Transcriptómica: solo tienen expresión en sangre CYP2D6 y CYP3A5, al comprobar en el Atlas se ve que es que los CYPs se expresan en el hígado y no en sangre 
-Epi: se deja anotado en CpGs_por_CYP_ventana100kb_con_betas


### 25/05/26
### Objetivo
Sacar los modelos para la relacion entre gene counts y metilacion 
### Transcripcion y epi: Modelo simple: Mvalues ~ gene_counts (para todos los CYPs)
## 1.Transcriptomica: pasar los nombres de las muestras de gene counts a formato epi 
Con el script_relacion_nombres y el excel relacion_nombres_muestras_epi_trancriptomica que tiene ambos formatos, se cambia el excel CYP_gene_counts_nombres_transcriptomica_238_individuos a los nombres de las muestras de epi CYP_gene_counts_nombres_epi_238_individuos
## 2.Epi: calcular los M-values a partir de los B-values para cada CYP
Con el script_transformar_Bvalues_a_Mvalues se pasan los Bvalues del excel CpGs_por_CYP_ventana100kb_con_betas a Mvalues (CpGs_por_CYP_ventana100kb_con_Mvalues.xlsx)
Se usa la formula M=log2(B/(1-B))
## 3.Fusion: Aplicar modelos 
Se adapta el script de Ines del pipeline de epigenomica (14_CaseControl_Analysis_Females_Males_limma_blood_prueba_sin_CtrlProbes3y4) a nuestro caso script modelos_epi_transcriptomica 
### Resultados 
-Se sacan cg significativas en algunos genes 

### 27/05/26
### Objetivo
Seguir con modelos_epi_transcriptomica
-Ver como poner el nombre para ver que cg es en el modelo 
-Correr el modelo simple para los demas CYPS 
-Analizar las significativas
-Ver si es necesario el eBayes o si serviaria solo el lmfit en nuestro caso
-Ir metiendo covariables (ver si se pueden usar los lambda para ver que modelo es el mejor)
### Transcripcion y epi: Modelo simple: Mvalues ~ gene_counts (para todos los CYPs)
Se analizan las cg significativas, 
1. Se considera significativas si tienen un P-valor-ajustado (FDR) >0.05 y solo salen tres para CYP3A5 
-Con LogFC positivo: cg02084114  (más expresión = más metilación)
-Con LogFC negativo: cg12694063 y cg13596235 (más expresión = menos metilación)
2. Se aplica el Bonferrroni (P-valor < 0.05/CpGs) y solo salen dos de los tres de antes de CYP3A5
-Con LogFC positivo: cg02084114  (más expresión = más metilación)
-Con LogFC negativo: cg12694063 (más expresión = menos metilación)
### Transcripcion y epi:  Modelo 1:  Mvalues ~ gene_counts + age + sex (para CYP2D6, CYP3A5  porque son las que tienen gene counts)
Sale lo mismo que el modelo simple
### Transcripcion y epi: Modelo 2:  Mvalues ~ gene_counts + age + sex + smoking (para CYP2D6 y CYP3A5 porque son las que tienen gene counts)
Sale lo mismo que el modelo simple
### Transcripcion y epi: Modelo 3:  Mvalues ~ gene_counts + age + sex + smoking + caso/control (para CYP2D6 y CYP3A5 porque son las que tienen gene counts)
Se analizan las cg significativas, 
1. Se considera significativas si tienen un P-valor-ajustado (FDR) >0.05 y solo sale una para CYP3A5 
-Con LogFC negativo: cg12694063 (más expresión = menos metilación)
2. Se aplica el Bonferrroni (P-valor < 0.05/CpGs) y sigue saliendo el de antes de CYP3A5
-Con LogFC negativo: cg12694063 (más expresión = menos metilación)
### Transcripcion y epi: Modelo 4:  Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 (para CYP2D6 y CYP3A5 porque son las que tienen gene counts)
No sale nada para ninguno de los genes 
### Transcripcion y epi: Modelo completo:  Mvalues ~ gene_counts + age + sex + smoking + caso/control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 (para CYP2D6 y CYP3A5 porque son las que tienen gene counts)
No sale nada para ninguno de los genes 
### Resulatados
Todo para CYP3A5 con Bonferroni
-MS: Mvalues ~ gene_counts: cg02084114, cg12694063
-M1: Mvalues ~ gene_counts + age + sex: cg02084114, cg12694063
-M2: Mvalues ~ gene_counts + age + sex + smoking: cg02084114, cg12694063
-M3: Mvalues ~ gene_counts + age + sex + smoking + caso/control: cg12694063
-M4: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: NADA 
-MC: Mvalues ~ gene_counts + age + sex + smoking + caso/control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: NADA 
De momento nos quedamos con el M2 (Mvalues ~ gene_counts + age + sex + smoking)


### 28/05/26
### Objetivo
-Aplicar sobre CYP3A5 M5, M6 y M7 de Transcripcion y epi
-Mirar que son las cg significativas de CYP3A5
### Resultados 
Todo para CYP3A5 con Bonferroni
-MS: Mvalues ~ gene_counts: cg02084114, cg12694063
-M1: Mvalues ~ gene_counts + age + sex: cg02084114, cg12694063
-M2: Mvalues ~ gene_counts + age + sex + smoking: cg02084114, cg12694063
-M3: Mvalues ~ gene_counts + age + sex + smoking + caso/control: cg12694063
-M4: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: NADA 
-M5: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2: NADA 
-M6: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2: cg02084114, cg12694063
-M7: Mvalues ~ gene_counts + age + sex + smoking + anc_1 + anc_2: cg02084114, cg12694063
-M8: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: cg02084114, cg12694063
-M9: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + caso/control: cg12694063
-M10: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + cell_1: NADA 
-MC: Mvalues ~ gene_counts + age + sex + smoking + caso/control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: NADA 
De momento nos quedamos con el M6 (Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2) porque tienen variables biologicas y tecnicas, se quita cell porque nos elimina los cg significativos, caso/control porque no es lo que buscamos con nuestro estudio de los CYPs y ancestriaporque la poblacion es toda europea y no es el objetivo de nuestro estudio y hemos visto que no afecta, asi que nos quedamos con lo minimo neceario para que sea fiable y no haya sobreajuste  


### 29/05/26
### Objetivo
-Hacer el modelo de nuevo sobre las gene_counts normalizadas  
-Mirar que son las cg significativas de CYP3A5
### Transcripcion y epi: Modelo simple: Mvalues ~ gene_counts (para todos los CYPs)
Solo salen significativas por Bonferroni en 
## Para CYP3A5
Se considera significativas si tienen un P-valor-ajustado (FDR) >0.05 y Bonferrroni (P-valor < 0.05/CpGs)
-Con LogFC positivo: cg02084114  (más expresión = más metilación)
-Con LogFC negativo: cg12694063, cg13596235, cg04706801 y cg05008948 (más expresión = menos metilación)
## Para CYP3A4
Se considera significativas si tienen un P-valor-ajustado (FDR) >0.05 y Bonferrroni (P-valor < 0.05/CpGs)
-Con LogFC positivo: cg11449766  (más expresión = más metilación)
-Con LogFC negativo: cg09914773 (más expresión = menos metilación)
## Para los demas CYPs no sale nada 
### Modelos con covariables
## Para CYP3A5 
-MS: Mvalues ~ gene_counts: cg02084114, cg12694063, cg13596235, cg04706801 y cg05008948 (5cg)
-M1: Mvalues ~ gene_counts + age + sex: cg02084114, cg12694063, cg13596235, cg04706801 y cg05008948 (5cg)
-M2: Mvalues ~ gene_counts + age + sex + smoking: cg02084114, cg12694063, cg13596235 y cg04706801 (4cg)
-M3: Mvalues ~ gene_counts + age + sex + smoking + caso/control: cg12694063 y cg13596235 (2cg)
-M4: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: cg12694063 (1cg)
-M5: Mvalues ~ gene_counts + age + sex + smoking + cell_1 + cell_2: cg12694063 (1cg)
-M6: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2: cg02084114, cg12694063, cg13596235 y cg04706801 (4cg)
-M7: Mvalues ~ gene_counts + age + sex + smoking + anc_1 + anc_2: cg02084114, cg12694063, cg13596235 y cg04706801 (4cg) 
-M8: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: cg02084114, cg12694063, cg13596235 y cg04706801 (4cg) 
-M9: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 + caso/control: cg12694063 y cg13596235 (2cg)
-MC: Mvalues ~ gene_counts + age + sex + smoking + caso/control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: cg12694063 (1cg)
De momento nos quedamos con el M8: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2 porque tienen variables biologicas y tecnicas, se quita cell porque nos elimina los cg significativos, caso/control porque no es lo que buscamos con nuestro estudio de los CYPs (se mete ancestria aunque hemos visto que no afecta porque son todo europeaos), asi que nos quedamos con lo minimo neceario para que sea fiable y no haya sobreajuste  
## Para CYP3A4 
-MS: Mvalues ~ gene_counts: cg11449766 y cg09914773 (2cg)
-M1: Mvalues ~ gene_counts + age + sex: cg11449766 y cg09914773 (2cg)
-M6: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2: cg11449766 y cg09914773 (2cg)
-MC: Mvalues ~ gene_counts + age + sex + smoking + caso/control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: cg11449766 y cg09914773 (2cg)
Da igual el modelo que se coja salen siempre esas dos CpGs asi que tambien M8 
### cg significativas de CYP3A5 
Se analizan las cg significativas en CpGs_por_CYP_ventana100kb_con_betas (para ver donde estan) y se marcan en amarillo 
1. Se considera significativas si tienen un P-valor-ajustado (FDR) >0.05 
-cg02084114 (LogFC positivo: más expresión = más metilación): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99337507%2D99937506&hgsid=4031321989_w4YZAy5JodQVbAyoa1WZezXWYnbx
-cg12694063 (LogFC negativo: más expresión = menos metilación): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99531818%2D99731817&hgsid=4031331559_IbEHtV4kM0g20nywdJV1gauUmiQQ
-cg13596235 (LogFC negativo: más expresión = menos metilación): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99532028%2D99732027&hgsid=4031502761_oK6msl8c1aKahG7mw89Ip6Hkia0d
-cg04706801 (LogFC negativo: más expresión = menos metilación): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99524259%2D99724258&hgsid=4031935157_CCFgVliD8gFAA4ReVQ4xDHn8wyZo
-cg05008948 (LogFC negativo: más expresión = menos metilación): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99533591%2D99733590&hgsid=4031937755_SlqfxAtdPHat5athFcJjbC6sYNGq (esta se va al aplicar el M8)
### cg significativas de CYP3A4 
Se analizan las cg significativas en CpGs_por_CYP_ventana100kb_con_betas (para ver donde estan) y se marcan en amarillo 
-cg11449766 (LogFC positivo: más expresión = más metilación): no sale buscando la cg en UCSC asi que se busca su posicion (99860727) cae en CYP3A43 https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99754539%2D99868513&hgsid=4031942255_kRKx7wyuAudr4gAQB4vW6a1371AY
-cg09914773 (LogFC negativo: más expresión = menos metilación): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99752211%2D99789608&hgsid=4031942255_kRKx7wyuAudr4gAQB4vW6a1371AY



### 01/06/26
### Objetivo
-Transcripcion y metabolizador: Hacer regresion lineal tipo metabolizador con expresion (gene counts) sin covariables (gene_counts ~ metabolizador (factor) SOLO para los CYPS con info de metabolizador)
### Transcripcion y metabolizador: Modelo simple sobre los 238 individuos de epi 
## 1.Transcriptomica: ya estan las gene_counts (CYP_gene_counts_normalizadas_nombres_transcriptomica_238_individuos) 
## 2.Tipo metabolizador: filtrar y formatear los datos de input 
Se parte del excel salida_genostar_filtrado_epi y se filtra para solo tener la info del matbolizer status y de la misma forma, se guarda en CYP_metabolizador_nombres_genostar_238_individuos
## 3.Fusion: Aplicar modelos 
Con script_modelos_metabolizador_transcriptomica, se hace un modelo lineal sin basarse en limma y no hay correlación en ninguno de los casos
No sale nada significativo 
-CYP2D6: 0.7955
-CYP2B6: 0.1235
-CYP3A5: 0.4165
-CYP2C9: 0.8315
-CYP2C19: 0.9055
### Resultados 
-Transcripcion y metabolizador:No sale nada significativo para regresion transcriptomica y tipo de metabolizador para los 238 individuos de epi 



### 02/06/26
### Objetivo
-Transcripcion y metabolizador: y Hacer regresion lineal tipo metabolizador con expresion (gene counts) sin covariables (gene_counts ~ metabolizador (factor) SOLO para los CYPS con info de metabolizador) para los 660 individuos 
-Transcripticion y epi: Mirar si las CpGs del CYP3A5 son DMR mirando si el P-valor (no ajustado) hace una escalera
-Trascripcion y epi: Hacer la regresion lineal para el ZSCAN5 en base a las cg significativas de CYP3A5
### Transcripcion y metabolizador: Modelo simple sobre los 660 individuos de transcriptomica 
## 1.Transcriptomica: filtrar las gene counts para los 660 individuos de genostar (CYP_gene_counts_normalizadas_nombres_transcriptomica_660_individuos)
Se parte de gene_counts_722_individuos_normalizados y se filtra para los 660 individuos de genostar y para los CYPS 
## 2.Tipo metabolizador: filtrar y formatear los datos de input 
Se parte del excel salida_genostar_filtrado_transcriptomica y se filtra para solo tener la info del matbolizer status y de la misma forma, se guarda en CYP_metabolizador_nombres_genostar_660_individuos
## 3.Fusion: Aplicar modelos 
Con script_modelos_metabolizador_transcriptomica, se hace un modelo lineal sin basarse en limma y no hay correlación en ninguno de los casos
No sale nada significativo 
-CYP2D6: 0.6132
-CYP2B6: 0.3147
-CYP3A5: 0.2816
-CYP2C9: 0.5997
-CYP2C19: 0.8067
### Transcripcion y epi CYP3A5: analizar posible DMR 
Analizando las posiciones de las 5 CpGs signifivativas se observa que cg02084114, cg12694063 y cg13596235 estan muy juntas y pueden ser una DMR y hay que comprobarlo 
## 1. Obetener los p-valores (no ajustados) de las CpGs de CYP3A5 
Se sacan con el script script_modelos_epi_transcriptomica y se guardan en DMR_CYP3A5
## 2. Obtener las posicion de las CpGs 
En en el excel DMR_CYP3A5 se hace BUSCARV a la hoja CYP3A5 del excel CpGs_por_CYP_ventana100kb_con_betas para obetener las posiciones 
## 3. Se analiza si las CpGs cercanas hacen una escalera en sus p-valores (Ballon plot)
En el excel DMR_CYP3A5 se ordena por posicion y se busca la zona con las CpGs (se marcan en rojo) y hace un ballon plot con la herramienta 
Link: http://www.bioinformatics.com.cn/plot_basic_ballon_plot_048_en
- Usuario: encuestadesatisfaccioncbmso@gmail.com
- Contraseña: 123456
### Transcripcion y epi: modelo simple: Mvalues ~ gene_counts (para ZSCAN5)
Se sigue el procedimiento del dia 25/05/26
## 1. Transcriptomica: gene counts del gen ZSCAN25
1. gene_counts_722_individuos_normalizados
2. ZSCAN25_gene_counts_normalizadas_nombres_transcriptomica_722_individuos
3. ZSCAN25_gene_counts_normalizadas_nombres_transcriptomica_238_individuos
4. ZSCAN25_gene_counts_normalizadas_nombres_epi_238_individuos (TRANSCIPRTOMICA)
## 2. Epi: 104 CpGs del CYP3A5
CpGs_por_CYP_ventana100kb_con_Mvalues.xlsx lo uncioque se va a usar solo la hoja de CYP3A5
## 3.Fusion: Aplicar modelos 
Script modelos_epi_transcriptomica al final en apartado ZSCAN25 y se saca el excel DMPs_ZSCAN25_MS 
## 4. Evaluar si las DMPs de ZSACN25 estan en la zona DMRs de CYP3A5 
En el excel DMPs_ZSCAN25_MS se hace BUSCARV a la hoja Zona DMR del excel DMRs_CYP3A5
### Resultados
-Transcripcion y metabolizador: No sale nada significativo para regresion transcriptomica y tipo de metabolizador para los 660 individuos de transcriptomica, tiene sentido porque el tipo de metabolizador no depende de la expresion del gen sino de la svariantes que tiene y esas son invariables, solo afectaria si tiene alguna variante framshift o de stop, que hace que no se produzca la enzima 
-DMR en CYP3A5 en el ballon plot 
-DMPs en ZSCAN25: queda por analizar 


### 03/06/26
-Poster CBM 


### 05/06/26
### Objetivo
-Ballon plot DRM en DMR_CYP3A5 (credenciales en el correo)
-Analizar DMPs en ZSCAN25 y su relacion con las de CYP3A5
### Transcripcion y epi CYP3A5: analizar posible DMR 
Analizando las posiciones de las 5 CpGs signifivativas se observa que cg02084114, cg12694063 y cg13596235 estan muy juntas y pueden ser una DMR y hay que comprobarlo 
## 3. Se analiza si las CpGs cercanas hacen una escalera en sus p-valores (Ballon plot)
En el excel DMR_CYP3A5 se ordena por posicion y se busca la zona con las CpGs (se marcan en rojo) y hace un ballon plot con la herramienta 
Link: http://www.bioinformatics.com.cn/plot_basic_ballon_plot_048_en
- Usuario: encuestadesatisfaccioncbmso@gmail.com
- Contraseña: 123456
Los p-valores transforman a -log10 para qu e el cambio sea proporcional (¿Por qué usamos el signo negativo en −log10(p‑value)?
Porque queremos que los p‑valores más pequeños (más significativos) se conviertan en números más grandes.
Si NO pusiéramos el signo negativo: Un p‑valor pequeño (0.0001) daría un número negativo grande (−4).  Un p‑valor grande (0.8) daría un número casi cero (−0.09).
Eso es lo contrario de lo que queremos en un gráfico donde el tamaño del globo representa significancia.)
Se guarda en Ballon_plot_DMR_CYP3A5
### Transcripcion y epi: modelo simple: Mvalues ~ gene_counts (para ZSCAN5)
## 4. Evaluar si las DMPs de ZSACN25 estan en la zona DMRs de CYP3A5 
En el excel DMPs_ZSCAN25_MS se hace BUSCARV a la hoja Zona DMR del excel DMRs_CYP3A5 y se hace hoja filtrada para FDR (quedan 29 cg) y para Bonferroni (quedan 13 cg). Sobre la hoja filtrada por Bonferroni se añade una columna con el p-valor de esas cg (BUSCARV a hoja DMR del excel DMRs_CYP3A5)
Resultado: De las cibnco significativas de CYP3A5, cuatro estan en el supuesto DMR y la otra cg05008948 no esta en DMR pero es significativo. Se analizan las CpGs de CYP3A5 para  el gen ZSACN25 y salen 13 significativas, de las cuales tres eran de las cinco signitivas de CYP3A5 (la suleta cg05008948 y dos del DMR cg02084114 y cg04706801), las 10 restantes no son significativas en CYP3A5 pero si para ZSACN25 y de esas diez, seis estan en el supuesto DMR de CYP3A5 y las otras cuatro no estan en el DMR de CYP3A5
### Transcripcion y epi: otros modelos para para ZSCAN5
## 1. Transcriptomica: gene counts del gen ZSCAN25
ZSCAN25_gene_counts_normalizadas_nombres_epi_238_individuos (TRANSCIPRTOMICA)
## 2. Epi: 104 CpGs del CYP3A5
CpGs_por_CYP_ventana100kb_con_Mvalues.xlsx lo uncioque se va a usar solo la hoja de CYP3A5
## 3.Fusion: Aplicar modelos para ZSCAN25
Script modelos_epi_transcriptomica al final en apartado ZSCAN25 y se saca el excel DMPs_ZSCAN25 
-MS: Mvalues ~ gene_counts: (13 cg) En DMPs_ZSCAN25_MS hoja filtrada por Bonferroni 
-M1: Mvalues ~ gene_counts + age + sex: 
-M6: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2: (14cg)
-M8: Mvalues ~ gene_counts + age + sex + smoking + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: (14cg) 
-MC: Mvalues ~ gene_counts + age + sex + smoking + caso/control + cell_1 + cell_2 + Ctrl_1 + Ctrl_2 + anc_1 + anc_2: 0cg
Nos quemos con el M8 porque asi se ha decidio antes 
## 4. Evaluar si las DMPs de ZSACN25 estan en la zona DMRs de CYP3A5 
En el excel DMPs_ZSCAN25_M8 se hace BUSCARV a la hoja Zona DMR del excel DMRs_CYP3A5 y se hace hoja filtrada para Bonferroni (quedan 14 cg). Sobre la hoja filtrada por Bonferroni se añade una columna con el p-valor de esas cg (BUSCARV a hoja DMR del excel DMRs_CYP3A5)
Resultado: De las cibnco significativas de CYP3A5, cuatro estan en el supuesto DMR y la otra cg05008948 no esta en DMR pero es significativo. Se analizan las CpGs de CYP3A5 para  el gen ZSACN25 y salen 14 significativas, de las cuales tres eran de las cinco signitivas de CYP3A5 (la suleta cg05008948 y dos del DMR cg02084114 y cg04706801), las 11 restantes no son significativas en CYP3A5 pero si para ZSACN25 y de esas once, seis estan en el supuesto DMR de CYP3A5 y las otras cinco no estan en el DMR de CYP3A5
### Resultados 
-DMR en CYP3A5 en Ballon_plot_DMR_CYP3A5
-DMPs en ZSCAN25: con el modelo simple hay 13 cgs que coinciden bastante con la supuesta DMR y con el M8 hay 14cgs que coinciden bastante con la supuesta DMR, por lo que se llega a la conclusion de que esa zona tiene una tendencia y puede ser una posible DMR  


### 08/06/26
### Objetivo
-Ver que metaboliza CYP3A5 
-Ver correo process CYP2D6, pero creo que el bajo call rate va a venir de que no tneemos rs que definen alelos segun proces_cyp2d6 y sobretodo se pierden muchos *1
### Genostar bajo call rate CYP2D6
Se ve que Genostar para asignar los diplotipos de CYP2D6 usa:
-como para los demas CYPs la tabla de referencia (Allele_def) que usa 150rs de los cuales solo tenemos 40rs, los 110 que nos faltan no son importantes
-por process_cyp2d6 que usa 15 SNPs, de los cuales solo tenemos 10rs, los 5 que faltan son importantes porque definen allel star especificos 
 -rs35742686: para *1 y *3. En gnomad frecuencia de 0.01750
 -rs5030655: para *1 y *6. En gnomad frecuencia de 0.01242
 -rs5030656: para *1 y *9. En gnomad frecuencia de 0.03065
 -rs5030865A: para *1, *8 y *14. En gnomad frecuencia de 0.00001531
 -rs5030865T: *1, *8 y *14. En gnomad frecuencia de 0.00001020
Se adapta la funciones process_cyp2d6_diana_adaptada para eliminar esos rs y se corren los scripts 
-script_3_genostar: se genera salida_genostar_filtrado_CYP2D6_adaptado.xlsx, CYP2D6_pie_chart_nuevo.png y CYP2D6_individuos_NA_nuevo.xlsx
-script_CALL_RATE: 97,12%
Pero al mirar las frecuencias de los metabolizadores no cuadra mucho con lo de Clinpgx por lo que el error que se esta cometiendo es alto aunque mejore mucho el Call rate 
### CYP3A5: fármacos metabolizados 
https://pmc.ncbi.nlm.nih.gov/articles/PMC11625447/
As indispensable members of the cytochrome P450 enzyme family, CYP3A4/5 are integral to the metabolism of a vast array of substrates, which can be categorized as follows: Steroid hormones, such as cortisol (Van Keulen et al., 2020), testosterone (Wu et al., 2022), and estradiol (Cheng & Zhou, 2000). Additionally, numerous lipophilic drugs necessitate CYP3A4/5 for metabolism, including benzodiazepines (Rodriguez-Antona et al., 2022), calcium channel blockers (Fuhr et al., 2022), antibiotics (Lenard et al., 2023), and specific antiepileptic medications (Zhao et al., 2021a), all of which are vital in the treatment of a myriad of diseases. While CYP3A4 and CYP3A5 share a broad substrate specificity, there are notable differences in their affinity and capacity to metabolize certain drugs. For instance, CYP3A4 is the primary enzyme involved in the metabolism of the immunosuppressant cyclosporine, which is crucial in post-organ transplant therapy (Zhai et al., 2024). In the metabolism of the immunosuppressant tacrolimus, where CYP3A5 is more efficient than CYP3A4 (Hannachi et al., 2024). 
Para las drugs que tiene recomendacion: https://www.clinpgx.org/gene/PA131/labelAnnotation



### 08/06/26
### Objetivo
-Hablar con Ines del mensaje de Claudio de las CpGs 
-Error CYP2D6 
### Error CYP2D6
Se deja apuntado en el cuaderno
### Boxplot DMR CYP3A5 
## 1. Epi: coger los Mvalues solo de las 4 CpGs significativas para CYP3A5 que estan en el DMR 
Se filtra CpGs_por_CYP_ventana100kb_con_Mvalues a CpGsig_DMR_CYP3A5_Mvalues
## 2. Transcriptomica: se coge CYP_gene_counts_normalizadas_nombres_epi_238_individuos 
## 3. Boxplots  
-Boxplot_DMR_CYP3A5_all_CpGs: con el script script_boxplot_DMR_CYP3A5_all_CpGs, que evalua la relacion entre el nivel de metilacion medio (se hace la media antes de todas las CpGs y luego ya la mediana para hacer el nivel) y las gene counts. No nos da mucha informacion porque junta las cuatro CpGs 
-Boxplot_DMR_CYP3A5_CpGs_separadas_mediana: con el script script_boxplot_DMR_CYP3A5_CpGs_separadas, evalua el nivel de metilacion (usando la mediana) y las gene counts para cada CpG significativa. Se ve que la CpG con LogFC positivo a más metilacion más gene counts y para los otros tres CpGs con LogFC negativo, a más metilación menos gene counts (se hace el Wilcoxon para ver significacia y todas son significativas excepto cg12694063)
-Boxplot_DMR_CYP3A5_CpGs_separadas_deciles: con el script script_boxplot_DMR_CYP3A5_CpGs_separadas, evalua el nivel de metilacion (usando los deciles (bajo para deciles 1,2,3 y 4 y alto para deciles 7,8,9 y 10)) y las gene counts para cada CpG significativa. Se ve que la CpG con LogFC positivo a más metilacion más gene counts y para los otros tres CpGs con LogFC negativo, a más metilación menos gene counts (se hace el Wilcoxon para ver significacia y todas son significativas)
### Resultados 
Los bvoxplots reafirman lo encontrado antes respecto a la DMR y el wilcoxon sale significativo 



-Evaluar los boxplots y losp-valores  
-Decidir el modelo para aplicar en transcriptomica y epi
-Hablar del errorde CYP2D6 con Ines

### Error CYP2D6 
El SNP rs1058164, incluido en algunas tablas de definición de alelos de CYP2D6 como modulador de expresión, no estaba disponible en nuestra plataforma y no se consideró en la asignación de fenotipos.


BIBLIOGRAFIA
https://pmc.ncbi.nlm.nih.gov/articles/PMC11625447/
As indispensable members of the cytochrome P450 enzyme family, CYP3A4/5 are integral to the metabolism of a vast array of substrates, which can be categorized as follows: Steroid hormones, such as cortisol (Van Keulen et al., 2020), testosterone (Wu et al., 2022), and estradiol (Cheng & Zhou, 2000). Additionally, numerous lipophilic drugs necessitate CYP3A4/5 for metabolism, including benzodiazepines (Rodriguez-Antona et al., 2022), calcium channel blockers (Fuhr et al., 2022), antibiotics (Lenard et al., 2023), and specific antiepileptic medications (Zhao et al., 2021a), all of which are vital in the treatment of a myriad of diseases. While CYP3A4 and CYP3A5 share a broad substrate specificity, there are notable differences in their affinity and capacity to metabolize certain drugs. For instance, CYP3A4 is the primary enzyme involved in the metabolism of the immunosuppressant cyclosporine, which is crucial in post-organ transplant therapy (Zhai et al., 2024). In the metabolism of the immunosuppressant tacrolimus, where CYP3A5 is more efficient than CYP3A4 (Hannachi et al., 2024). 
https://pmc.ncbi.nlm.nih.gov/articles/PMC10197210/pdf/dmd.122.001007.pdf

