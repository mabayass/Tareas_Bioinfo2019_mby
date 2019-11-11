
# GITHUB
## *Unidad 2 “Organización de un proyecto Bioinformático”*

***
***
GitHub es una plataforma informática que permite desarrollar proyectos de códigos fuente mediante un software utilizando el sistema de
control de versiones Git, siendo una de sus mayores ventajas la capacidad de colaboración de diferentes usuarios sobre un mismo proyecto 
a distancia. 

Durante las primeras semanas del segundo módulo en el curso avanzado de Tópicos computacionales en Genética y Genómica, se ha trabajado
en la Unidad 2 y 7 del workshop en github asociado a la materia AliciaMstt/BioinfinvRepro. La Unidad 2 da inicio con un tutorial general
sobre Markdown, cuyas aplicaciones se insertan en este avance del módulo II. La Unidad 2 permitió la creación de distintos repositorios con el fin de aplicar los conocimientos básicos de los tutoriales, donde encontramos el  [Repositorio 1](https://github.com/mabayass/Tareas_Bioinfo2019_mby) y el [Repositorio 2](https://github.com/mabayass/mirepointro) como los repositorios personales del práctico. A continuación, se observan imágenes del terminal de GitHub donde se aplicaron los comandos otorgados por los tutoriales como _Hello-world_ y _An Intro to Git and GitHub for begginers_.  

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/tutorial%20git%20imagen%201.png "Tutorial GitHub")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/tutorial%20git%20imagen%202.png)


***

***


# Expresión Diferencial
## *Unidad 7 “Análisis de Expresión diferencial a partir de secuencias de RNA”*

Se desarrolla un análisis bioinformático de expresión diferencial a partir de muestras de *Sulfolobus acidocaldarius* con sacuenciación de RNA, en el cual se compararon muestras mutadas y no mutadas en distintos medios de cultivo, lo que permitió determinar genes diferencialmente expresados producto de los diferentes estímulos al que fueron sometidos durante la experimentación, donde el análisis bioinformático cumple un papel importante en la interpretación de estos datos mediante gráficos y tablas.   

Para fines de este práctico, se utilziaron datos de 4 librerías de lecturas de arqueobacterias *Sulfolobus acidocaldarius*, organismo con un único cromosoma circular con 2,225,959 pares de bas, con un contenido de G+C de 36.7%. El gen asociado a la formación de biopelículas fue sometido a un *knockdown* para luego ser expuesto a diferentes medios de cultivo y asi poder estudiar los cambios en la expresión génica.

Las muestras fueron identificadas de la siguiente forma:

* Wild type P
* Wild type B
* Mutant P
* Mutant B

A continuación se enumeran los pasos determinados por el tutorial de la Unidad 7 para realizar el análisis de la expresion génica. 

### Pasos tutorial Unidad 7

#### 1. Crear variables
*Carpetas preexistentes que contienen la ubicación de las carpetas que ya se encuentran creadas dentro del directorio home*

RAW=/shared/bioinfo1/common/raw_data/

ANN=/shared/bioinfo1/common/annot/

REF=/shared/bioinfo1/common/ref_genome/

#### 2. Crear carpetas de salida 
*Carpetas con rutas donde se almacenará información, y que permiten guardar textos producto de los análisis bioinformáticos*

QC=../qc

FIL=../filtered

ALN=../alignment

CNT=../count

#### 3. Control de Calidad 
*Directorios donde se almacena información procesada. Luego se ejecuta el programa IlluQC_PRLL.pl, el cual genera un reporte completo de la calidad de las secuencias. La version PRLL permite ejecutar el comando usando distintos CPU al mismo tiempo. Una vez finalizado el análisis, el programa arroja distintos gráficos que representan la calidad de las lecturas, el contenido GC y otros datos necesarios para el análisis posterior de la secuenciación.* 

* Control de Calidad Muestra **Wild Type P**

Comando ejecutado en Unix
> illuqc -se "$RAW/MW001_P.fastq" 5 A -onlystat -t 2 -o "wild_planctonic" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_QualRangePerBase.png "Cantidad de lecturas por base")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_avgQual.png "Valor promedio de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_gcDistribution.png "Distribucion de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_qualDistribution.png "Distribución de calidad")


##### CONCLUSIÓN 


***


* Control de Calidad Muestra **Wild Type B**

Comando ejecutado en Unix
> illuqc -se "$RAW/MW001_B3.fastq" 5 A -onlystat -t 2 -o "wild_biofilm" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_QualRangePerBase.png "Cantidad de lecturas por base")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_gcDistribution.png "Distribución de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_qualDistribution.png "Distribución de calidad")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_summary.png "Resumen de calidad de lecturas")

#####CONCLUSIÓN

***


* Control de Calidad Muestra **Mutant P**

Comando ejecutado en Unix
> illuqc -se "$RAW/0446_P.fastq" 5 A -onlystat -t 2 -o "mut_planctonic" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_avgQual.png "Valor promedio de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_gcDistribution.png "Distribución de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_qualDistribution.png "Distribución de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_summary.png "Resumen de calidad de lecturas")

##### CONCLUSIÓN

***


* Control de Calidad Muestra **Mutant B**

Comando ejecutado en Unix
> illuqc -se "$RAW/0446_B3.fastq" 5 A -onlystat -t 2 -o "mut_biofilm" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_QualRangePerBase.png "Cantidad de lecturas por base")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_gcDistribution.png "Distribución de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_qualDistribution.png "Distribución de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_summary.png "Resumen de calidad de lecturas")

##### CONCLUSIÓN



***

#### 4. Filtro de secuencias
*Luego de obtener los resultados del control de calidad de la secuenciación de RNA, las librerías son filtradas con el objetivo de eliminar lecturas con calidad menor de 20% en el 80% de la extensión, cuyos resultados genera librerías de lectura que seran utilizadas en el Alineamiento de las secuencias.*

> Se crea un nuevo directorio _FIL_ con aquellas carpetas donde se almacenarán los resultados del proceso de filtrado. 

* Filtro de secuencias de **Wild Type P**
> illuqc -se "$RAW/MW001_P.fastq" 5 A -l 80 -s 20 -t 2 -o "wild_planctonic" -c 1 & 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/MW001_P.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***


* Filtro de secuencias de **Wild Type B**
>illuqc -se "$RAW/MW001_B3.fastq" 5 A -l 80 -s 20 -t 2 -o "wild_biofilm" -c 1

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/MW001_B3.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Mutant P**
>illuqc -se "$RAW/0446_P.fastq" 5 A -l 80 -s 20 -t 2 -o "mut_planctonic" -c 1

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/0446_P.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***


* Filtro de secuencias de **Mutant B**
>illuqc -se "$RAW/0446_B3.fastq" 5 A -l 80 -s 20 -t 2 -o "mut_biofilm" -c 1 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/0446_B3.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")


##### CONCLUSIÓN


***


#### 5. Alineamiento
*A partir de las librerias de lectura producidas a partir de la filtración del punto anterior, se procede a hacer un alineamiento de la secuenciación de RNA de las muestras frente al genoma de referencia. A continuación se muestran los comandos que permitieron el análisis.*

* Muestra Wild Type P 
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/wild_planctonic/MW001_P.fastq_filtered" > "ALN/MW001_P_aligned.sam" &

* Muestra Wild Type B
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/wild_biofilm/MW001_B3.fastq_filtered" > "ALN/MW001_B3_aligned.sam" &

* Muestra Mutant P
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/mut_planctonic/0446_P.fastq_filtered" > "ALN/0446_P_aligned.sam" &

* Muestra Wild Type B
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/mut_biofilm/0446_B3.fastq_filtered" > "ALN/0446_B3_aligned.sam" &

A partir de los comandos ejecutados, es posible observar los archivos tipo _.sam_ en la carpeta _ALN_, a continuación: 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ALN/alineados.png "Screenshot archivos .sam creados")


##### Característica de un archivo .sam 


***






#### 6. Estimación de la abundancia
*Para poder hacer una estimación de las lecturas mapeadas en cada uno de los genes, se debe instalar un programa llamado HTSeq-Count versión 0.6.1, cuyos archivos emitidos serán utilizados para el análisis de expresión diferencial.*

+ Instalación de HTSeq-Count versión 0.6.1 en el directorio code (de acuerdo al tutorial)
> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ **pip install HTseq

+ Uso del comando para estimación de lecturas mapeadas
> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ **python -m HTSeq.scripts.count -t Gene -i GenID "ALN/MW001_P_aligned.sam" "$ANN/saci.gff3" > "CNT/MW001_P.count" &

> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ **python -m HTSeq.scripts.count -t Gene -i GenID "ALN/MW001_B3_aligned.sam" "$ANN/saci.gff3" > "CNT/MW001_B3.count" &**

> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ **python -m HTSeq.scripts.count -t Gene -i GenID "ALN/0446_P_aligned.sam" "$ANN/saci.gff3" > "CNT/0446_P.count" &**

> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ **python -m HTSeq.scripts.count -t Gene -i GenID "ALN/0446_B3_aligned.sam" "$ANN/saci.gff3" > "CNT/0446_B3.count" &**


A partir de los comandos ejecutados, es posible observar los archivos tipo _.count_ en la carpeta _CNT_ a continuación:
![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/CNT/count.png "Screenshot archivos .count creados")


***








#### 7. Análisis de Expresión Diferencial
*El análisis de expresion diferencial es realizado en Rstudio, en el cual se pudo cargar el paquete edgeR exitosamente*

###### Preparación de los datos para el análisis de expresión diferencial
+ a. Crear directorios para almacenar gráficos y tablas de análisis

input_dir  <- **file.path("/Users/mjgr1/Documents/Master Genetica/Segundo Semestre/Bioinformatica/Unidad II","ExpDif")**
output_pseudo <- file.path("..","diff_expr", "pseudocounts")
output_histogram <- file.path("..","diff_expr", "histograms")
output_pvalue_fdr <- file.path("..","diff_expr", "pvalue_fdr")
output_table <- file.path("..","diff_expr", "tables")


+ b. Se crean las carpetas de salida luego de comprobar que las carpetas anteriormente creadas existen

if(!file.exists(input_dir)){stop("Data directory doesn't exist: ", input_dir)}
if(!file.exists(output_pseudo)){dir.create(output_pseudo, mode = "0755", recursive=T)}
if(!file.exists(output_histogram)){dir.create(output_histogram, mode = "0755", recursive=T)}
if(!file.exists(output_pvalue_fdr)){dir.create(output_pvalue_fdr, mode = "0755", recursive=T)}
if(!file.exists(output_table)){dir.create(output_table, mode = "0755", recursive=T)}


+ c. Cargar la librería 'edgeR' 
  
**library(edgeR)


+ d. Lectura de archivos _count_ obtenidos de Unix

**Wild type planctonic**
wild_p <- read.delim(file=file.path(input_dir, "MW001_P.count"), sep="\t", header = F, check=F); dim(wild_p);colnames(wild_p) <- c("Gen_ID", "Count")

**Wild type biofilm**
wild_b <- read.delim(file=file.path(input_dir, "MW001_B3.count"), sep="\t", header = F, check=F); dim(wild_b);colnames(wild_b) <- c("Gen_ID", "Count")

**Mutant planctonic**
mut_p <- read.delim(file=file.path(input_dir, "0446_P.count"), sep="\t", header = F, check=F); dim(mut_p); colnames(mut_p) <- c("Gen_ID", "Count")

**Mutant biofilm**
mut_b <- read.delim(file=file.path(input_dir, "0446_B3.count"), sep="\t", header = F, check=F); dim(mut_b); colnames(mut_b) <- c("Gen_ID", "Count")


+ e. Colapsar los datasets
rawcounts2 <- data.frame(wild_p$Gen_ID, WildType_P = wild_p$Count, WildType_B = wild_b$Count, Mutant_B = mut_b$Count, row.names = 1)


+ f. Remover las columnas que no seran usadas en el analisis
to_remove <- rownames(rawcounts2) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")


+ g. RPKM (Reads per kilobase per million o Lecturas por kilobase por millon)
rpkm <- cpm(rawcounts2)


+ h. Establecer las lecturas que serán utilizadas en el análisis y hacer un filtrado
keep <- rowSums(rpkm > 1) >= 3 & !to_remove
rawcounts2 <- rawcounts2[keep,]



###### Análisis de expresión diferencial por medios de cultivo

+ a. Crear vector para muestras agrupadas
group_culture <- c("planctonic","biofilm","biofilm")


+ b. Crear un objeto DGE (Expresión génica diferencial)
dge_culture <- DGEList(counts = rawcounts2, group = group_culture)


+ c. Normalizar factores por tamaño de librería
dge_culture <- calcNormFactors(dge_culture)


+ d. Estimar dispersion por muestra y por gen
dge_culture <- estimateCommonDisp(dge_culture)
dge_culture <- estimateTagwiseDisp(dge_culture)


+ e. Hacer el test Exact, cuyo análisis se basa en asumir conteos de distribución binomial negativa
de_culture <- exactTest(dge_culture, pair = c("planctonic","biofilm"))


+ f. Obtener resumen de los resultados
results_culture <- topTags(de_culture, n = nrow(dge_culture))
results_culture <- results_culture$table

+ g. Obtener ID de genes expresados diferencialmente por medio de cultivo
ids_culture <- rownames(results_culture[results_culture$FDR < 0.1,])




###### Análisis de expresión diferencial por genotipos
Se siguen los pasos que en el punto anterior

+ a. Crear un set de data de COUNTS sin genes expresados diferencialmente por medio de cultivo
rawcounts_genotype <- rawcounts2[!rownames(rawcounts2) %in% ids_culture,]


+ b. Crear un vector por muestras agrupadas
group_genotype <- c("wildtype","wildtype","mutant")


+ c. Crear un objeto DGE
dge_genotype <- DGEList(counts = rawcounts_genotype, group = group_genotype)


+ d. Normalizar factores por tama;o de libreria
dge_genotype <- calcNormFactors(dge_genotype)


+ e. Estimar dispersion por muestra y por gen
dge_genotype <- estimateCommonDisp(dge_genotype)
dge_genotype <- estimateTagwiseDisp(dge_genotype)


+ f. Hacer el test Exact 
de_genotype <- exactTest(dge_genotype, pair = c("wildtype","mutant"))


+ g. Obtener resumen de los resultados
results_genotype <- topTags(de_genotype, n = nrow(de_genotype))
results_genotype <- results_genotype$table


+ h. Obtener genes expresados diferencialmente por medio de cultivo
ids_genotype <- rownames(results_genotype)
ids_genotype <- ids_genotype[results_genotype$FDR < .1]



###### Generación de resultados

+ a. Establecer vectores Booleans para definir genes con expresion diferencial para ambos factores

1. Medio de cultivo
de_genes_culture <- rownames(rawcounts2) %in% ids_culture

2. Genotipo
de_genes_genotype <- rownames(rawcounts2) %in% ids_genotype


+ b. Obtener los pseudocounts obtenidos del exact test y transformarlos en escala logarítmica
pseudocounts <- data.frame(rownames(rawcounts2), WildType_P = log10(dge_culture$pseudo.counts[,1]), WildType_B = log10(dge_culture$pseudo.counts[,2]), Mutant_B = log10(dge_culture$pseudo.counts[,3]), DE_C = de_genes_culture, DE_G = de_genes_genotype, row.names = 1)

En este paso se resaltan los genes diferencialmente expresados 


+ c. Gráficos y archivos PDF con expresión diferencial según **medios de cultivo

pdf(file=file.path(output_pseudo,"pair_expression_culture.pdf"), width = 8, height = 4)
par(mfrow = c(1,2))

plot(pseudocounts$WildType_P, pseudocounts$WildType_B, col = ifelse(pseudocounts$DE_C, "red", "blue"), main = "Wild Type", xlab = "Planctonic", ylab = "Biofilm", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
abline(lsfit(pseudocounts$WildType_P, pseudocounts$WildType_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/WT%20ambos%20medios.png "Gráfico WT Medios de cultivo")


+ d. Gráficos y archivos PDF con expresión diferencial según **Genotipo

pdf(file=file.path(output_pseudo,"pair_expression_genotype.pdf"), width = 8, height = 4)
par(mfrow = c(1,2))

plot(pseudocounts$WildType_P, pseudocounts$Mutant_B, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "Planctonic", xlab = "Wild Type", ylab = "Mutant", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
abline(lsfit(pseudocounts$WildType_P, pseudocounts$Mutant_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/con%20abline.png "Gráfico WT-Mut Planctonic")

plot(pseudocounts$WildType_B, pseudocounts$Mutant_B, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "Biofilm", xlab = "Wild Type", ylab = "Mutant", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
abline(lsfit(pseudocounts$WildType_B, pseudocounts$Mutant_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/WTB%20con%20MUTb.png "Gráfico WT-Mut Biofilm")

plot(pseudocounts$WildType_P, pseudocounts$WildType_B, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "WildType", xlab = "Wild Type Planctonic", ylab = " WT BIOFILM", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
abline(lsfit(pseudocounts$WildType_P, pseudocounts$WildType_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/entre%20wt.png "Gráfico Wild Type")


+ e. Creación de Histogramas de los P-values

pdf(file=file.path(output_histogram,"histograms_pvalue.pdf"), width = 8, height = 4)
par(mfrow = c(1,2))

hist(x = results_culture$PValue, col = "skyblue", border = "blue", main = "Culture", xlab = "P-value", ylab = "Frequency", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/HISTOGRAMA%20MEDIOS.png "Histograma de Medios de Cultivo")

hist(x = results_genotype$PValue, col = "skyblue", border = "blue", main = "Genotype", xlab = "P-value", ylab = "Frequency", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/histograma%20genotipo.png "Histograma de Genotipo")


+ f. Gráfico de P-value vs FDR
pdf(file=file.path(output_pvalue_fdr, "pvalue_fdr.pdf"), width = 8, height = 4)
par(mfrow = c(1,2))

plot(results_culture$PValue, results_culture$FDR, col = "blue", main = "Culture", xlab = "P-value", ylab = "FDR", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/FDR%20medios.png "Plot FDR vs P-Value de Medios de Cultivo")

plot(results_genotype$PValue, results_genotype$FDR, col = "blue", main = "Genotype", xlab = "P-value", ylab = "FDR", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/FDR%20genotipos.png "Plot FDR vs P-Value de Genotipo")


+ g. Resumen en Tabla de resultados

1. Medio de cultivo
write.table(x=results_culture, file=file.path(output_table, "table_de_genes_culture.csv"), quote=F, sep="\t", dec=".", row.names=T, col.names=T)

2. Genotipo
write.table(x=results_genotype, file=file.path(output_table, "table_de_genes_genotype.csv"), quote=F, sep="\t", dec=".", row.names=T, col.names=T)

