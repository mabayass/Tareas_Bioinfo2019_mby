
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
> illuqc -se "GM1598-1_R1_001.fastq" 5 A -onlystat -t 2 -o "GM1598-1_R1_001" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-1_R1_001.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-1_R1_001.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-1_R1_001.fastq_gcDistribution.png "Distribucion de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-1_R1_001.fastq_qualDistribution.png "Distribución de calidad")


***


* Control de Calidad Muestra **Input Infectado**

Comando ejecutado en Unix
> illuqc -se "GM1598-2_R1_001.fastq" 5 A -onlystat -t 2 -o "GM1598-2_R1_001" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-2_R1_001.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-2_R1_001.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-2_R1_001.fastq_gcDistribution.png "Distribucion de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-2_R1_001.fastq_qualDistribution.png "Distribución de calidad")


***


* Control de Calidad Muestra **Pulldown Basal**

Comando ejecutado en Unix
> illuqc -se "GM1598-5_R1_001.fastq" 5 A -onlystat -t 2 -o "GM1598-5_R1_001" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-5_R1_001.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-5_R1_001.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-5_R1_001.fastq_gcDistribution.png "Distribucion de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-5_R1_001.fastq_qualDistribution.png "Distribución de calidad")


***


* Control de Calidad Muestra **Pulldown Infectado**

Comando ejecutado en Unix
> illuqc -se "GM1598-6_R1_001.fastq" 5 A -onlystat -t 2 -o "GM1598-6_R1_001" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-6_R1_001.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-6_R1_001.fastq_baseCompostion.png "Composición de Bases")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-6_R1_001.fastq_gcDistribution.png "Distribucion de Contenido de GC")

***

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-6_R1_001.fastq_qualDistribution.png "Distribución de calidad")


***

##### CONCLUSIÓN

Todas las muestras bajo las distintas condiciones presentan una secuenciación que cumple con los parámetros de calidad. Para el caso del Quality score range vemos que se encuentran dentro del Q30, el cual significa que la probabilidad de que la identificación de casa base secuenciada presente un error de 1 en 1000, el cual cumple con lo establecido para la secuenciación con HiSeq. Tanto el porcentaje de las bases nucleotídicas como la distribucion de CG cumple con lo esperado para muestras de RNA en _archeas_. 


***
***


### 4. Filtro de secuencias
*Luego de obtener los resultados del control de calidad de la secuenciación de RNA, las librerías son filtradas con el objetivo de eliminar lecturas con calidad menor de 20 (Q20) en el 80% de la extensión, cuyos resultados genera librerías de lectura que seran utilizadas en el Alineamiento de las secuencias.*


* Filtro de secuencias de **Input Basal**
> illuqc -se "GM1598-1_R1_001.fastq" 5 A -l 80 -s 20 -t 2 -o "GM1598-1_R1_001" -c 1 & 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-1_R1_001.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Input Infectado**
> illuqc -se "GM1598-2_R1_001.fastq" 5 A -l 80 -s 20 -t 2 -o "GM1598-2_R1_001" -c 1 & 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-2_R1_001.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Pulldown Basal**
> illuqc -se "GM1598-5_R1_001.fastq" 5 A -l 80 -s 20 -t 2 -o "GM1598-5_R1_001" -c 1 & 1

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-5_R1_001.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")

***

* Filtro de secuencias de **Pulldown Infectado**
> illuqc -se "GM1598-6_R1_001.fastq" 5 A -l 80 -s 20 -t 2 -o "GM1598-6_R1_001" -c 1 & 

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/GM1598-6_R1_001.fastq_filtered_QualRangePerBase.png "Cantidad de lecturas por base")


##### CONCLUSIÓN

Se obtiene una filtración exitosa que genera archivos que permiten alinear las muestras con el genoma de referencia. 


***
***


### 5. Alineamiento
*A partir de las librerias de lectura producidas a partir de la filtración del punto anterior, se procede a hacer un alineamiento de la secuenciación de RNA de las muestras frente al genoma de referencia. A continuación se muestran los comandos que permitieron el análisis. En este caso, se utiliza el comando bwa (Burrow Wheelers Aligment) dado que permite tomar el genoma de referencia y alinearlo correctamente con las muestras en estudio. Sin embargo, existen otros comandos como bowtie2 que ejecuta la misma acción generando archivos tipo .sam para el paso siguiente que sería la estimación de la abundancia.*

*Los archivos .sam consiste en archivos de texto delimitados con tabs con información general tanto de la secuenciación como del alineamiento que presenta 11 campos obligatorios y otros opcionales, donde cada renglón de la sección de alineamientos corresponde a un segmento de lectura de cada fragmento introducido en el secuenciador.*


Los comandos llevados a cabo para el alineamiento:

* Muestra **Input Basal** 
> bwa078 mem "$REF/genome.fasta" -t 1 "GM1598-1_R1_001.fastq_filtered" > "GM1598-1_aligned.sam" &

* Muestra **Input Infectado**
> bwa078 mem "$REF/genome.fasta" -t 1 "GM1598-2_R1_001.fastq_filtered" > "GM1598-2_aligned.sam" &

* Muestra **Pulldown Basal**
> bwa078 mem "$REF/genome.fasta" -t 1 "GM1598-5_R1_001.fastq_filtered" > "GM1598-5_aligned.sam" &

* Muestra **Pulldown Infectado**
> bwa078 mem "$REF/genome.fasta" -t 1 "GM1598-6_R1_001.fastq_filtered" > "GM1598-6_aligned.sam" &

Los archivos _.sam_ generados son los archivos de entrada para la estimación de la abundancia mediante el programa HTSeq que se muestra a continuación:


***
***



### 6. Estimación de la abundancia
*Para poder hacer una estimación de las lecturas mapeadas en cada uno de los genes, se debe instalar un programa llamado HTSeq-Count versión 0.6.1, cuyos archivos emitidos serán utilizados para el análisis de expresión diferencial. Este análisis se realiza a partir del archivo .sam emitido en el paso anterior donde una vez ejecutado el comando, te arroja un archivo .count que permite hacer posteriormente el análisis de expresión diferencial, siendo la estimación de la abundancia el último paso de importancia durante el análisis de control de calidad de la secuenciación.*


+ Uso del comando para estimación de lecturas mapeadas

* Muestra **Input Basal** 
> /usr/bin/python -m HTSeq.scripts.count -t Gene -i gene_id "GM1598-1_aligned.sam" "Homo_sapiens.GRCh38.96.chr" > "GM1598-1.count" &

* Muestra **Input Infectado**
> /usr/bin/python -m HTSeq.scripts.count -t Gene -i gene_id "GM1598-2_aligned.sam" "Homo_sapiens.GRCh38.96.chr" > "GM1598-2.count" &

* Muestra **Pulldown Basal**
> /usr/bin/python -m HTSeq.scripts.count -t Gene -i gene_id "GM1598-5_aligned.sam" "Homo_sapiens.GRCh38.96.chr" > "GM1598-5.count" &

* Muestra **Pulldown Infectado**
> /usr/bin/python -m HTSeq.scripts.count -t Gene -i gene_id "GM1598-6_aligned.sam" "Homo_sapiens.GRCh38.96.chr" > "GM1598-6.count" &



A partir de los comandos ejecutados se obtienen los archivos requeridos para poder llevar a cabo el análisis de expresión diferencial.


***
***



### 7. Análisis de Expresión Diferencial
El análisis de DE se lleva a cabo en el programa Rstudio

#### Preparación de los datos para el análisis de expresión diferencial
+ a. Crear directorios para almacenar gráficos y tablas de análisis

> input_dir  <- **file.path("/Users/mjgr1/Documents/Master Genetica/Segundo Semestre/Bioinformatica/Unidad III","Count")**
> output_pseudo <- file.path("..","Count", "pseudocounts")
> output_histogram <- file.path("..","Count", "histograms")
> output_pvalue_fdr <- file.path("..","Count", "pvalue_fdr")
> output_table <- file.path("..","Count", "tables")

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
> library(ggplot2)


+ d. Lectura de archivos _count_ obtenidos de Linux en el paso del Control de Calidad.

**Wild type planctonic**
> inputbasal <- read.delim(file=file.path(input_dir, "GM1598-1.count"), sep="\t", header = F, check=F); dim(wild_p);colnames(inputbasal) <- c("Gen_ID", "Count")

**Wild type biofilm**
> inputinfectado <- read.delim(file=file.path(input_dir, "GM1598-2.count"), sep="\t", header = F, check=F); dim(wild_b);colnames(inputinfectado) <- c("Gen_ID", "Count")

**Mutant planctonic**
> pdbasal <- read.delim(file=file.path(input_dir, "GM1598-5.count"), sep="\t", header = F, check=F); dim(mut_p); colnames(pdbasal) <- c("Gen_ID", "Count")

**Mutant biofilm**
> pdinfectado <- read.delim(file=file.path(input_dir, "GM1598-6.count"), sep="\t", header = F, check=F); dim(mut_b); colnames(pdinfectado) <- c("Gen_ID", "Count")



+ e. Colapsar los datasets
> rawcounts <- data.frame(inputbasal$Gen_ID, InputBasal = inputbasal$Count, InputInfectado = inputinfectado$Count, PulldownBasal = pdbasal$Count, PulldownInfectado = pdinfectado$Count, row.names = 1)


+ f. Remover las columnas que no seran usadas en el analisis
> to_remove <- rownames(rawcounts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")


+ g. RPKM (Reads per kilobase per million o Lecturas por kilobase por millon)
> rpkm <- cpm(rawcounts)


+ h. Establecer las lecturas que serán utilizadas en el análisis y hacer un filtrado
> keep <- rowSums(rpkm > 1) >= 3 & !to_remove
> rawcounts <- rawcounts[keep,]


***


### Análisis de expresión diferencial por Condicion

+ a. Crear vector para muestras agrupadas
> group_condicion <- c("Basal","Basal","Infectado","Infectado")


+ b. Crear un objeto DGE (Expresión génica diferencial)
> dge_condicion <- DGEList(counts = rawcounts, group = group_condicion)


+ c. Normalizar factores por tamaño de librería
> dge_condicion <- calcNormFactors(dge_condicion)


+ d. Estimar dispersion por muestra y por gen
> dge_condicion <- estimateCommonDisp(dge_condicion)
> dge_condicion <- estimateTagwiseDisp(dge_condicion)


+ e. Hacer el test Exact, cuyo análisis se basa en asumir conteos de distribución binomial negativa
> de_condicion <- exactTest(dge_condicion, pair = c("Basal","Infectado"))


+ f. Obtener resumen de los resultados
> results_condicion <- topTags(de_condicion, n = nrow(dge_condicion))
> results_condicion <- results_condicion$table


+ g. Obtener ID de genes expresados diferencialmente por condicion
> ids_condicion <- rownames(results_condicion[results_condicion$FDR < 0.1,])


En este paso se lleva a cabo el análisis de expresión diferencial pero es posteriormente que se ejecutan los comandos para hacer gráficos y tablas para poder visualizar los datos

***



### Análisis de expresión diferencial por input/pulldown
Se siguen los pasos que en el punto anterior

+ a. Crear vector para muestras agrupadas
> group_entrada <- c("input","input","pulldown","pulldown")


+ b. Crear un objeto DGE (Expresión génica diferencial)
> dge_entrada <- DGEList(counts = rawcounts, group = group_entrada)


+ c. Normalizar factores por tamaño de librería
> dge_entrada <- calcNormFactors(dge_entrada)


+ d. Estimar dispersion por muestra y por gen
> dge_entrada <- estimateCommonDisp(dge_entrada)
> dge_entrada <- estimateTagwiseDisp(dge_entrada)


+ e. Hacer el test Exact, cuyo análisis se basa en asumir conteos de distribución binomial negativa
> de_entrada <- exactTest(dge_entrada, pair = c("input","pulldown"))


+ f. Obtener resumen de los resultados
> results_entrada <- topTags(de_entrada, n = nrow(dge_entrada))
> results_entrada <- results_entrada$table


+ g. Obtener ID de genes expresados diferencialmente por condicion
> ids_entrada <- rownames(results_entrada[results_entrada$FDR < 0.1,])


+ h. Obtener genes expresados diferencialmente por entrada
> ids_entrada <- rownames(results_entrada)
> ids_entrada <- ids_entrada[results_entrada$FDR < .1]


***


### Generación de resultados

+ a. Establecer vectores Booleans para definir genes con expresion diferencial para ambos factores

**1. Condicion**
> de_genes_condicion <- rownames(rawcounts) %in% ids_condicion

**2. Entrada**
> de_genes_entrada <- rownames(rawcounts) %in% ids_entrada



+ b. Obtener los pseudocounts obtenidos del exact test y transformarlos en escala logarítmica
> pseudocounts <- data.frame(rownames(rawcounts), InputBasal = log10(dge_condicion$pseudo.counts[,1]), InputInfectado = log10(dge_condicion$pseudo.counts[,2]), PulldownBasal = log10(dge_condicion$pseudo.counts[,3]), PulldownInfectado = log10(dge_condicion$pseudo.counts[,4]), DE_C = de_genes_condicion, DE_G = de_genes_entrada, row.names = 1)

_En este paso se resaltan los genes diferencialmente expresados_ 


+ c. Gráficos y archivos PDF con expresión diferencial según **Condicion

_Comando para creación de un documento PDF con gráfico de expresión diferencial entre condiciones_

> pdf(file=file.path(output_pseudo,"pair_expression_condicion.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comando para la creación de un gráfico de expresión diferencial entre las dos condiciones de input_

> plot(pseudocounts$InputBasal, pseudocounts$InputInfectado, col = ifelse(pseudocounts$DE_C, "red", "blue"), main = "Input", xlab = "Infectado", ylab = "Basal", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$InputBasal, pseudocounts$InputInfectado), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/WT%20ambos%20medios.png "Gráfico WT Medios de cultivo")

El plot muestra en puntos rojos aquellos genes que se encuentran diferencialmente expresados, lo que nos indica que bajo distintas condiciones de medios de cultivo, se activan o se reprimen genes asociados a la respuesta ambiental. 




+ d. Gráficos y archivos PDF con expresión diferencial según **Entrada

_Comando para creación de un documento PDF con gráfico de expresión diferencial entre entradas_

> pdf(file=file.path(output_pseudo,"pair_expression_entrada.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comando para la creación de un gráfico de expresión diferencial entre las dos entradas del análisis_

> plot(pseudocounts$PulldownBasal, pseudocounts$PulldownInfectado, col = ifelse(pseudocounts$DE_C, "red", "blue"), main = "Entrada", xlab = "Input", ylab = "Pulldown", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$PulldownBasal, pseudocounts$PulldownInfectado), col = "black")

> pdf(file=file.path(output_pseudo,"pair_expression_entrada.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))



_Comandos para la creación de gráficos de expresión diferencial entre las condiciones y las entradas infectado

> plot(pseudocounts$InputInfectado, pseudocounts$PulldownInfectado, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "Infectado", xlab = "Input", ylab = "Pulldown", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$InputInfectado, pseudocounts$PulldownInfectado), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/con%20abline.png "Gráfico WT-Mut Planctonic")


_Comandos para la creación de gráficos de expresión diferencial entre las condiciones y las entradas input

> plot(pseudocounts$InputBasal, pseudocounts$PulldownBasal, col = ifelse(pseudocounts$DE_G, "red", "blue"), main = "Basal", xlab = "Input", ylab = "Pulldown", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)
> abline(lsfit(pseudocounts$InputBasal, pseudocounts$PulldownBasal), col = "black")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/WTB%20con%20MUTb.png "Gráfico WT-Mut Biofilm")




+ e. Creación de Histogramas de los P-values

_Comando para creación de un documento PDF con histograma de P-value_

> pdf(file=file.path(output_histogram,"histograms_pvalue.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comandos para la creación de gráficos de histogramas de los genotipos y los medios de cultivo_

**1. Condicion** 
> hist(x = results_condicion$PValue, col = "skyblue", border = "blue", main = "Condicion", xlab = "P-value", ylab = "Frequency", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/HISTOGRAMA%20MEDIOS.png "Histograma de Medios de Cultivo")


**2. Entrada**
> hist(x = results_entrada$PValue, col = "skyblue", border = "blue", main = "Entrada", xlab = "P-value", ylab = "Frequency", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/histograma%20genotipo.png "Histograma de Genotipo")



+ f. Gráfico de P-value vs FDR

_Comando para creación de un documento PDF con gráfico FDR vs P-value_

> pdf(file=file.path(output_pvalue_fdr, "pvalue_fdr.pdf"), width = 8, height = 4)
> par(mfrow = c(1,2))


_Comandos para la creación de gráficos FDR vs P-value_

**1. Medios de cultivo**
> plot(results_condicion$PValue, results_condicion$FDR, col = "blue", main = "Condicion", xlab = "P-value", ylab = "FDR", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/FDR%20medios.png "Plot FDR vs P-Value de Medios de Cultivo")


**2. Entrada**
> plot(results_entrada$PValue, results_entrada$FDR, col = "blue", main = "Entrada", xlab = "P-value", ylab = "FDR", cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.2, las = 01)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/DE/FDR%20genotipos.png "Plot FDR vs P-Value de Genotipo")


El estadístico FDR (False Discovery Rate) determina aquellos falsos positivos dentro de los p-values significativos, por ende al graficar ambos análisis estamos modificando el p-value y restringiendo su valor para obtener un rechazo aceptación de a hipótesis nula mucho mas robusta. 


+ g. Resumen en Tabla de resultados

**1. Condicion**
> write.table(x=results_condicion, file=file.path(output_table, "table_de_genes_condicion.csv"), quote=F, sep="\t", dec=".", row.names=T, col.names=T)

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/tabla%20medios.png "Imagen Excel Medio de cultivo")


**2. Entrada**
> write.table(x=results_entrada, file=file.path(output_table, "table_de_genes_entrada.csv"), quote=F, sep="\t", dec=".", row.names=T, col.names=T)


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/tabla%20genotipo.png "Imagen Excel Genotipo")


***
***

##Conclusión 


