
***
***


# Expresión Diferencial
## *Unidad 7 “Análisis de Expresión diferencial a partir de secuencias de RNA”*

Bajo la importancia que presenta mantener el foco de atención en los estudios del VIH-1 por la incidencia que permanece en nuestra sociedad, se busca profundizar los conocimientos del tema asociando los lncARN y su efecto sobre la replicacion del virus. En este ensayo se pretende llevar a cabo un análisis bioinformático que determine los cambios de expresión de los lncARN bajo dos condiciones: con infección de VIH-1 y en condiciones basales, llevadas a cabo en células J-Lat, es decir, células Jurkat que fueron transformadas con un vector antirretroviral con un corrimiento del marco de lectura y una adición de GFP. Estas células permanecen en latencia hasta su activación con PMA, siendo este un activador de linfocitos T. Ambas condiciones fueron sometidas a un crooslinking con formaldehido y luego incubadas con sondas biotiniladas específicas de ARN viral, y luego purificadas con perlas de estreptavidina (que presentan afinidad por la biotina de las sondas). El extracto purificado se somete a elucion y Proteinasa K para luego secuenciar los ARN remanentes. Las muestras fueron secuenciadas por HiSeq con una lectura de 100 pb paired end realizadas en el servicio de secuenciación de la Universidad Mayor. 

Las muestras fueron identificadas de la siguiente forma:

* Input Basal (no infectado)
* Input Infectado
* Pulldown Basal (no infectado)
* Pulldown infectado

A continuacion se enumeran los pasos para realizar el análisis de expresión diferencial:

### 1. Control de Calidad 
*Directorios donde se almacena información procesada. Luego se ejecuta el programa IlluQC_PRLL.pl, el cual genera un reporte completo de la calidad de las secuencias. La version PRLL permite ejecutar el comando usando distintos CPU al mismo tiempo. Una vez finalizado el análisis, el programa arroja distintos gráficos que representan la calidad de las lecturas, el contenido GC y otros datos necesarios para el análisis posterior de la secuenciación.* 

* Control de Calidad Muestra **Input Basal**

Comando ejecutado en Unix
> illuqc -se "GM1598-1_R1_001.fastq" 5 A -onlystat -t 2 -o "GM1598-1_R1" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_avgQual.png "Valor promedio de calidad")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_gcDistribution.png "Distribucion de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_qualDistribution.png "Distribución de calidad")


***


* Control de Calidad Muestra **Wild Type B**

Comando ejecutado en Unix
> illuqc -se "$RAW/MW001_B3.fastq" 5 A -onlystat -t 2 -o "wild_biofilm" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_gcDistribution.png "Distribución de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_qualDistribution.png "Distribución de calidad")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_summary.png "Resumen de calidad de lecturas")


***


* Control de Calidad Muestra **Mutant P**

Comando ejecutado en Unix
> illuqc -se "$RAW/0446_P.fastq" 5 A -onlystat -t 2 -o "mut_planctonic" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_avgQual.png "Valor promedio de calidad")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_gcDistribution.png "Distribución de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_qualDistribution.png "Distribución de calidad")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_summary.png "Resumen de calidad de lecturas")


***


* Control de Calidad Muestra **Mutant B**

Comando ejecutado en Unix
> illuqc -se "$RAW/0446_B3.fastq" 5 A -onlystat -t 2 -o "mut_biofilm" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_gcDistribution.png "Distribución de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_qualDistribution.png "Distribución de calidad")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_summary.png "Resumen de calidad de lecturas")

***

##### CONCLUSIÓN

Todas las muestras genotípicas bajo las distintas condiciones presentan una secuenciación que cumple con los parámetros de calidad. Para el caso del Quality score range vemos que se encuentran dentro del Q30, el cual significa que la probabilidad de que la identificación de casa base secuenciada presente un error de 1 en 1000, el cual cumple con lo establecido para la secuenciación con HiSeq 2500. Tanto el porcentaje de las bases nucleotídicas como la distribucion de CG cumple con lo esperado para muestras de RNA en _archeas_. 


***
***


### 4. Filtro de secuencias
*Luego de obtener los resultados del control de calidad de la secuenciación de RNA, las librerías son filtradas con el objetivo de eliminar lecturas con calidad menor de 20 (Q20) en el 80% de la extensión, cuyos resultados genera librerías de lectura que seran utilizadas en el Alineamiento de las secuencias.*

Se crea un nuevo directorio _FIL_ con aquellas carpetas donde se almacenarán los resultados del proceso de filtrado. 


* Filtro de secuencias de **Wild Type P**
> illuqc -se "$RAW/MW001_P.fastq" 5 A -l 80 -s 20 -t 2 -o "wild_planctonic" -c 1 & 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/MW001_P.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Wild Type B**
> illuqc -se "$RAW/MW001_B3.fastq" 5 A -l 80 -s 20 -t 2 -o "wild_biofilm" -c 1

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/MW001_B3.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Mutant P**
> illuqc -se "$RAW/0446_P.fastq" 5 A -l 80 -s 20 -t 2 -o "mut_planctonic" -c 1

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/0446_P.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Mutant B**
> illuqc -se "$RAW/0446_B3.fastq" 5 A -l 80 -s 20 -t 2 -o "mut_biofilm" -c 1 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/FIL/0446_B3.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")


##### CONCLUSIÓN

Se obtiene una filtración exitosa que genera archivos que permiten alineas las muestras con el genoma de referencia. 


***
***


### 5. Alineamiento
*A partir de las librerias de lectura producidas a partir de la filtración del punto anterior, se procede a hacer un alineamiento de la secuenciación de RNA de las muestras frente al genoma de referencia. A continuación se muestran los comandos que permitieron el análisis. En este caso, se utiliza el comando bwa (Burrow Wheelers Aligment) dado que permite tomar el genoma de referencia y alinearlo correctamente con las muestras en estudio. Sin embargo, existen otros comandos como bowtie2 que ejecuta la misma acción generando archivos tipo .sam para el paso siguiente que sería la estimación de la abundancia.*

*Los archivos .sam consiste en archivos de texto delimitados con tabs con información general tanto de la secuenciación como del alineamiento que presenta 11 campos obligatorios y otros opcionales, donde cada renglón de la sección de alineamientos corresponde a un segmento de lectura de cada fragmento introducido en el secuenciador.*


Los comandos llevados a cabo para el alineamiento:

* Muestra **Wild Type P** 
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/wild_planctonic/MW001_P.fastq_filtered" > "ALN/MW001_P_aligned.sam" &

* Muestra **Wild Type B**
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/wild_biofilm/MW001_B3.fastq_filtered" > "ALN/MW001_B3_aligned.sam" &

* Muestra **Mutant P**
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/mut_planctonic/0446_P.fastq_filtered" > "ALN/0446_P_aligned.sam" &

* Muestra **Mutant B**
> bwa078 mem "$REF/genome.fasta" -t 1 "QC/FIL/mut_biofilm/0446_B3.fastq_filtered" > "ALN/0446_B3_aligned.sam" &



A partir de los comandos ejecutados, es posible observar los archivos tipo _.sam_ en la carpeta _ALN_, a continuación: 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ALN/alineados.png "Screenshot archivos .sam creados")


***
***



### 6. Estimación de la abundancia
*Para poder hacer una estimación de las lecturas mapeadas en cada uno de los genes, se debe instalar un programa llamado HTSeq-Count versión 0.6.1, cuyos archivos emitidos serán utilizados para el análisis de expresión diferencial. Este análisis se realiza a partir del archivo .sam emitido en el paso anterior donde una vez ejecutado el comando, te arroja un archivo .count que permite hacer posteriormente el análisis de expresión diferencial, siendo la estimación de la abundancia el último paso de importancia durante el análisis de control de calidad de la secuenciación.*


+ Instalación de HTSeq-Count versión 0.6.1 en el directorio code (de acuerdo al tutorial)
> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ pip install HTseq

+ Uso del comando para estimación de lecturas mapeadas
> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ python -m HTSeq.scripts.count -t Gene -i GenID "ALN/MW001_P_aligned.sam" "$ANN/saci.gff3" > "CNT/MW001_P.count" &

> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ python -m HTSeq.scripts.count -t Gene -i GenID "ALN/MW001_B3_aligned.sam" "$ANN/saci.gff3" > "CNT/MW001_B3.count" &

> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ python -m HTSeq.scripts.count -t Gene -i GenID "ALN/0446_P_aligned.sam" "$ANN/saci.gff3" > "CNT/0446_P.count" &

> bioinfo1@genoma1:~/mbarreto/TareaU7/code$ python -m HTSeq.scripts.count -t Gene -i GenID "ALN/0446_B3_aligned.sam" "$ANN/saci.gff3" > "CNT/0446_B3.count" &



A partir de los comandos ejecutados, es posible observar los archivos tipo _.count_ en la carpeta _CNT_ a continuación:
![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/CNT/count.png "Screenshot archivos .count creados")


***
***



### 7. Análisis de Expresión Diferencial
*El análisis de expresion diferencial es realizado en Rstudio, en el cual se pudo cargar el paquete edgeR exitosamente*

#### Preparación de los datos para el análisis de expresión diferencial
+ a. Crear directorios para almacenar gráficos y tablas de análisis

> input_dir  <- **file.path("/Users/mjgr1/Documents/Master Genetica/Segundo Semestre/Bioinformatica/Unidad II","ExpDif")**
> output_pseudo <- file.path("..","diff_expr", "pseudocounts")
> output_histogram <- file.path("..","diff_expr", "histograms")
> output_pvalue_fdr <- file.path("..","diff_expr", "pvalue_fdr")
> output_table <- file.path("..","diff_expr", "tables")

Con los comandos mencionados se establece la carpeta donde estan los archivos count emitidos en la estimacion de la abundancia y se crean directorios para crear carpetas de salida donde serán almacenados los gráficos del análisis. 

+ b. Se crean las carpetas de salida luego de comprobar que las carpetas anteriormente creadas existen

> if(!file.exists(input_dir)){stop("Data directory doesn't exist: ", input_dir)}
> if(!file.exists(output_pseudo)){dir.create(output_pseudo, mode = "0755", recursive=T)}
> if(!file.exists(output_histogram)){dir.create(output_histogram, mode = "0755", recursive=T)}
> if(!file.exists(output_pvalue_fdr)){dir.create(output_pvalue_fdr, mode = "0755", recursive=T)}
> if(!file.exists(output_table)){dir.create(output_table, mode = "0755", recursive=T)}

Con los comandos if! se ejecuta una verificación de creación de directorios y se crean entonces las carpetas de salida. 

+ c. Cargar la librería 'edgeR' 
> library(edgeR)


+ d. Lectura de archivos _count_ obtenidos de Linux en el paso del Control de Calidad.

**Wild type planctonic**
> wild_p <- read.delim(file=file.path(input_dir, "MW001_P.count"), sep="\t", header = F, check=F); dim(wild_p);colnames(wild_p) <- c("Gen_ID", "Count")

**Wild type biofilm**
> wild_b <- read.delim(file=file.path(input_dir, "MW001_B3.count"), sep="\t", header = F, check=F); dim(wild_b);colnames(wild_b) <- c("Gen_ID", "Count")

**Mutant planctonic**
> mut_p <- read.delim(file=file.path(input_dir, "0446_P.count"), sep="\t", header = F, check=F); dim(mut_p); colnames(mut_p) <- c("Gen_ID", "Count")

**Mutant biofilm**
> mut_b <- read.delim(file=file.path(input_dir, "0446_B3.count"), sep="\t", header = F, check=F); dim(mut_b); colnames(mut_b) <- c("Gen_ID", "Count")



+ e. Colapsar los datasets
> rawcounts2 <- data.frame(wild_p$Gen_ID, WildType_P = wild_p$Count, WildType_B = wild_b$Count, Mutant_B = mut_b$Count, row.names = 1)


+ f. Remover las columnas que no seran usadas en el analisis
> to_remove <- rownames(rawcounts2) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")


+ g. RPKM (Reads per kilobase per million o Lecturas por kilobase por millon)
> rpkm <- cpm(rawcounts2)


+ h. Establecer las lecturas que serán utilizadas en el análisis y hacer un filtrado
> keep <- rowSums(rpkm > 1) >= 3 & !to_remove
> rawcounts2 <- rawcounts2[keep,]


***


### Análisis de expresión diferencial por medios de cultivo

+ a. Crear vector para muestras agrupadas
> group_culture <- c("planctonic","biofilm","biofilm")


+ b. Crear un objeto DGE (Expresión génica diferencial)
> dge_culture <- DGEList(counts = rawcounts2, group = group_culture)


+ c. Normalizar factores por tamaño de librería
> dge_culture <- calcNormFactors(dge_culture)


+ d. Estimar dispersion por muestra y por gen
> dge_culture <- estimateCommonDisp(dge_culture)
> dge_culture <- estimateTagwiseDisp(dge_culture)


+ e. Hacer el test Exact, cuyo análisis se basa en asumir conteos de distribución binomial negativa
> de_culture <- exactTest(dge_culture, pair = c("planctonic","biofilm"))


+ f. Obtener resumen de los resultados
> results_culture <- topTags(de_culture, n = nrow(dge_culture))
> results_culture <- results_culture$table


+ g. Obtener ID de genes expresados diferencialmente por medio de cultivo
> ids_culture <- rownames(results_culture[results_culture$FDR < 0.1,])


En este paso se lleva a cabo el análisis de expresión diferencial pero es posteriormente que se ejecutan los comandos para hacer gráficos y tablas para poder visualizar los datos

***



### Análisis de expresión diferencial por genotipos
Se siguen los pasos que en el punto anterior

+ a. Crear un set de data de COUNTS sin genes expresados diferencialmente por medio de cultivo
> rawcounts_genotype <- rawcounts2[!rownames(rawcounts2) %in% ids_culture,]


+ b. Crear un vector por muestras agrupadas
> group_genotype <- c("wildtype","wildtype","mutant")


+ c. Crear un objeto DGE
> dge_genotype <- DGEList(counts = rawcounts_genotype, group = group_genotype)


+ d. Normalizar factores por tama;o de libreria
> dge_genotype <- calcNormFactors(dge_genotype)


+ e. Estimar dispersion por muestra y por gen
> dge_genotype <- estimateCommonDisp(dge_genotype)
> dge_genotype <- estimateTagwiseDisp(dge_genotype)


+ f. Hacer el test Exact 
> de_genotype <- exactTest(dge_genotype, pair = c("wildtype","mutant"))


+ g. Obtener resumen de los resultados
> results_genotype <- topTags(de_genotype, n = nrow(de_genotype))
> results_genotype <- results_genotype$table


+ h. Obtener genes expresados diferencialmente por medio de cultivo
> ids_genotype <- rownames(results_genotype)
> ids_genotype <- ids_genotype[results_genotype$FDR < .1]


***


### Generación de resultados

+ a. Establecer vectores Booleans para definir genes con expresion diferencial para ambos factores

**1. Medio de cultivo**
> de_genes_culture <- rownames(rawcounts2) %in% ids_culture

**2. Genotipo**
> de_genes_genotype <- rownames(rawcounts2) %in% ids_genotype



+ b. Obtener los pseudocounts obtenidos del exact test y transformarlos en escala logarítmica
> pseudocounts <- data.frame(rownames(rawcounts2), WildType_P = log10(dge_culture$pseudo.counts[,1]), WildType_B = log10(dge_culture$pseudo.counts[,2]), Mutant_B = log10(dge_culture$pseudo.counts[,3]), DE_C = de_genes_culture, DE_G = de_genes_genotype, row.names = 1)

_En este paso se resaltan los genes diferencialmente expresados_ 



+ c. Gráficos y archivos PDF con expresión diferencial según **medios de cultivo

_Comando para creación de un documento PDF con gráfico de expresión diferencial entre medios de cultivo_

> pdf(file=file.path(output_pseudo,"pair_expression_culture.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))



_Comando para la creación de un gráfico de expresión diferencial entre las dos condiciones de medios de cultivo en las muestras Wild Type_

> plot(pseudocounts$WildType_P, pseudocounts$WildType_B, col = ifelse(pseudocounts$DE_C, "red", "blue"), main = "Wild Type", xlab = "Planctonic", ylab = "Biofilm", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$WildType_P, pseudocounts$WildType_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/WT%20ambos%20medios.png "Gráfico WT Medios de cultivo")

El plot muestra en puntos rojos aquellos genes que se encuentran diferencialmente expresados, lo que nos indica que bajo distintas condiciones de medios de cultivo, se activan o se reprimen genes asociados a la respuesta ambiental. 




+ d. Gráficos y archivos PDF con expresión diferencial según **Genotipo

_Comando para creación de un documento PDF con gráfico de expresión diferencial entre genotipos_

> pdf(file=file.path(output_pseudo,"pair_expression_genotype.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comandos para la creación de gráficos de expresión diferencial entre los genotipos y los medios de cultivo_

> plot(pseudocounts$WildType_P, pseudocounts$Mutant_B, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "Planctonic", xlab = "Wild Type", ylab = "Mutant", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$WildType_P, pseudocounts$Mutant_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/con%20abline.png "Gráfico WT-Mut Planctonic")


> plot(pseudocounts$WildType_B, pseudocounts$Mutant_B, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "Biofilm", xlab = "Wild Type", ylab = "Mutant", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$WildType_B, pseudocounts$Mutant_B), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/WTB%20con%20MUTb.png "Gráfico WT-Mut Biofilm")

En los gráficos anteriores se observa ambos genotipos (wild type y mutante o knockout) a las diferentes condiciones experimentales, donde se observan algunos genes que se expresan en ambos genotipos bajo la misma condición, lo que explica que esos genes no se encuentran regulados por el gen Lrs14-like, el cual es el gen en estudio. 





+ e. Creación de Histogramas de los P-values

_Comando para creación de un documento PDF con histograma de P-value_

> pdf(file=file.path(output_histogram,"histograms_pvalue.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comandos para la creación de gráficos de histogramas de los genotipos y los medios de cultivo_

**1. Medios de cultivo** 
> hist(x = results_culture$PValue, col = "skyblue", border = "blue", main = "Culture", xlab = "P-value", ylab = "Frequency", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/HISTOGRAMA%20MEDIOS.png "Histograma de Medios de Cultivo")


**2. Genotipos**
> hist(x = results_genotype$PValue, col = "skyblue", border = "blue", main = "Genotype", xlab = "P-value", ylab = "Frequency", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/histograma%20genotipo.png "Histograma de Genotipo")


En el histograma de medios de cultivo se observa una mayor frecuencia de genes con un p-value significativo, lo que confirma los resultados obtenidos en los plots anteriores. Por el contrario, el histograma de genotipos no presenta una tendencia sino mas bien, genes que en su mayoría, se encuentran con p-values estadísticamente no significativos. 



+ f. Gráfico de P-value vs FDR

_Comando para creación de un documento PDF con gráfico FDR vs P-value_

> pdf(file=file.path(output_pvalue_fdr, "pvalue_fdr.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comandos para la creación de gráficos FDR vs P-value_

**1. Medios de cultivo**
> plot(results_culture$PValue, results_culture$FDR, col = "blue", main = "Culture", xlab = "P-value", ylab = "FDR", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/FDR%20medios.png "Plot FDR vs P-Value de Medios de Cultivo")


**2. Genotipo**
> plot(results_genotype$PValue, results_genotype$FDR, col = "blue", main = "Genotype", xlab = "P-value", ylab = "FDR", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/FDR%20genotipos.png "Plot FDR vs P-Value de Genotipo")


El estadístico FDR (False Discovery Rate) determina aquellos falsos positivos dentro de los p-values significativos, por ende al graficar ambos análisis estamos modificando el p-value y restringiendo su valor para obtener un rechazo aceptación de a hipótesis nula mucho mas robusta. 


+ g. Resumen en Tabla de resultados

**1. Medio de cultivo**
> write.table(x=results_culture, file=file.path(output_table, "table_de_genes_culture.csv"), quote=F, sep="\t", dec=".", row.names=T, col.names=T)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/tabla%20medios.png "Imagen Excel Medio de cultivo")


**2. Genotipo**
> write.table(x=results_genotype, file=file.path(output_table, "table_de_genes_genotype.csv"), quote=F, sep="\t", dec=".", row.names=T, col.names=T)


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/tabla%20genotipo.png "Imagen Excel Genotipo")


***
***

##Conclusión 

El análisis de expresión diferencial permitió, en primer lugar, concluir que 165 genes presentaron una expresión distinta a la basal cuando se somete la archea _Sulfolobus acidocaldarius_ a diferentes condiciones ambientales, sea este un medio plantónico, donde las arqueas flotan libre en eun medio acuoso, o sea un medio biofilm donde gracias a un sustrato, estos microorganismos puedan adherirse y ensamblarse en una estructura que permite formar la biopelícula. Por otro lado, solo un gen, el gen Saci_1295, presentó una expresion diferencial cuando se compararon los genotipos en cada condicion. Esto implica que la mutacion en el gen Lrs14-like, pudiese estar regulando otras vías que no implica respuestas a cambios ambientales, dado que bajo la condición mutada no hubo respuestas diferenciales significativas. 

El proseguir del estudio pudiese direccionarse a realizar réplicas del experimento con el fin de reproducir los resultados y verificar que en efecto se obtienen resultados fieles o a definir por ejemplo, experimentos que impliquen aquellos genes diferenciados en los medios de cultivo para identificar con mayor precision aquelas vías de expresión génica implicada en el metabolismo de las archeas frente a estímulos ambientales. 
