# Proyecto Diana
## 19/02/26

### Objetivo
Obtener los bfile de los 660 individuos con las variantes SNPs y SNVs de los genes CYP (CYP1A2, CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP3A4 y CYP3A5)

### Pasos
#### 1. Sacar de los 5374 individuos, los 722 individuos con info de RNA 
Se eliminan 62 porque no estan en la tabla de los 5374 (ascendencia no europea u otras razones)

Resultado: 660indv.txt (familia (columna 1), ID individuo (columna 2))

#### 2. Sacar la info necesaria de los CYPs  
cromosoma, start y end de UCSC Genome Browser (a start y end se le pone una ventana de 5000)

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
```bash
awk '{print $4}' 

```


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
Se filtra de TODOS los individuos (df_CYP3A4.csv) solo los rs que aparecern en la tabla de referencia que usa genostar, sacada de github genostar (CYP3A5_Allele_def.rda). Todo esta en el script filtrar_rs_referencia_genostar.
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP3A4_Allele_def_rs_excluidos)
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP3A4_Allele_def_rs_excluidos)
#### Correr el genostar de TODOS: assign_diplotype 
Se saca el diplotipo de loS 660 individuos (CYP3A5_resultados_final_diplotipos)

#### Correr el genostar de TODOS: all_geno_pheno 
Da un error al intentar sacar el metabolizador 

### Resultados
Para CYP1A2 ya se tienen los diplotipos de 609 individuos en CYP1A2_resultados_final.rds y los 51 individuos restantes en CYP1A2_resultados_final_NA.rds, todos sacados con genostar.  Para el CYP3A4 ya se tienen los diplotipos de  620 individuos en CYP3A4_resultados_final.rds y de los 40 individuos restantes en CYP3A4_resultados_final_NA.rds, todos sacados con genostar. 
Para CYP1A2 ya se tienen los diplotipos de los 660 individuos en CYP3A5_resultados_final_diplotipos.rds, pero no los metabolizadores 



## 26/03/26
### Objetivo
Reunión con Claudio 
### Resultados
Me dijo que hasta el momento bien, pero que antes de seguir, me dijo de ir gen por gen y ver porque había rs que nosotros tenemos genotipados pero no usa genostar y ver porque genostar usa rs que no tenemos nosotros. Si es información redundante o si se esta dejando información importante fuera 


## 30/03/26
### Objetivo
Analizar rs de CYP1A2  
### Resumen rs CYP1A2
Hay 22 rs en comun entre genostar y los nuestros. Nosotros tenemos 60 más que genostar no usa. Genostar usa 24 más que nosotros no tenemos 
#### Tabla 
Se crae la tabla CYP1A2_TABLA_COMPLETA con la info de los rs 
#### LIMITACIONES: guardar los rs que usa genostar pero no estan en nuestros datos 
Todo esta en el script filtrar_rs_referencia_genostar. Esto se va a poder usar de limitación para ver si hay posibles alelos que no estemos viendo al hacer el filtrado (CYP3A4_Allele_def_rs_excluidos)

### Resultados
Me dijo que hasta el momento bien, pero que antes de seguir, me dijo de ir gen por gen y ver porque había rs que nosotros tenemos genotipados pero no usa genostar y ver porque genostar usa rs que no tenemos nosotros. Si es información redundante o si se esta dejando información importante fuera 
