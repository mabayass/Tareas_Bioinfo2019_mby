
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
> bwa078 mem "Homo_sapiens.GRCh38.96.chr" -t 1 "GM1598-1_R1_001.fastq_filtered" > "GM1598-1_aligned.sam" &

* Muestra **Input Infectado**
> bwa078 mem "Homo_sapiens.GRCh38.96.chr" -t 1 "GM1598-2_R1_001.fastq_filtered" > "GM1598-2_aligned.sam" &

* Muestra **Pulldown Basal**
> bwa078 mem "Homo_sapiens.GRCh38.96.chr" -t 1 "GM1598-5_R1_001.fastq_filtered" > "GM1598-5_aligned.sam" &

* Muestra **Pulldown Infectado**
> bwa078 mem "Homo_sapiens.GRCh38.96.chr" -t 1 "GM1598-6_R1_001.fastq_filtered" > "GM1598-6_aligned.sam" &

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
El análisis de DE se lleva a cabo en el programa R.

#### Preparación de los datos para el análisis de expresión diferencial
+ a. cat counts.txt | Rscript deseq1.r 

+ b. Se leen los argumentos de la linea de comandos
> args = commandArgs(trailing0nly=TRUE)
if (lenght(args)!=1) {stop("Experimental design must be specified as: NxM at the command line", call= FALSE)}
first = args[1]
if (first == 'install') {source("http://bioconductor.org/biocLite.R")
   biocLite("DESeq")
   stop("Installation completed", call.+FALSE)}

+ c. Extraer el dise;o experimental de la linea de comando
design = unlist(strplit(first, 'x'))

+ d. Encontrar los counts del diseño
cond1_num = as.integer(desing[1])
cond2_num = as.integer(desing[2])

+ e. Establecer las condiciones basadas en el setup experimental
cond_1= rep("cond1", cond1_num)
cond_2 = rep("cond2", cond2_num)

+ f. Cargar la librería
library(DESeq)

+ g. Cargar la data del archivo
counts= read.table("stdin", header=TRUE, row.names=1. sep="\t")

+ h. Dado que Kallisto genera unos counts estaimados como numeros reales y DESeq1 permite solo numeros enteros se necesita convertir numeros reales en enteros
int_counts = as.matrix(counts)
int_counts = apply(int_counts, 2, function(x) as.integer(round(x)))
row.names(int_counts) <-row.names(counts)

+ i. Reemplazar con counts enteros
counts <- int_counts

+ j. Establecer condiciones
conditions = factor(c(cond_1, cond_2))

+ k. Crear una tabla de los counts
cds = newCountDataSet(counts, conditions)

+ l. Estimar tamaño de los factores
cds = estimateSizeFactors(cds)

+ m. Estimar dispersiones
cds = estimateDispersions(cds)

+ n. Comparar 
results = nbiomTest(cds, "cond1", "cond2")

+ o. Hacer un sort de la tabla de resultados de las columnas de padj y foldChange
sorted = results[with(results, order(padj, -foldChange)),]

+ p. Escribir los resultados del standard output
write.table(sorted, file="", sep="\t", row.name=FALSE, quote=FALSE)

+ q. Mantener solo los valores diferencialmente expresados
diffs <- subset(sorted, padj < 0.05, select=c(id))

+ r. Normalizar los counts y escribir en un archivo
nc = counts(cds, normalized=TRUE)

+ s. Transformarlo en un dataframe para gener nombres de columnas adecuadas
dt = data.frame("id"=rownames(nc), nc)

+ t. Mantener solo las filas que presentaban expresion diferencial antes
keep = subset(dt, id %in% diffs$id)

+ u. Guardar en la matrix de data normalizada
write.table(keep, file="norm-matrix-deseq1.txt", sep="\t", row.name=FALSE, col.names=TRUE, quote=FALSE) 



***


### Gráficos

+ a. Contar las caracteristicas de cada nombre de  los genes
> featureCounts --primary -C -t transcript -p -T 12 -g transcript_id -a Homo_sapiens.GRCh38.96.chr.gtf -o countspulldown.txt    pulldowninfeccion_n1.bam pulldowninfeccion_n2.bam pulldownmock_n1.bam pulldownmock_n2.bam  2>> logpulldown.txt
featureCounts --primary -C -t transcript -p -T 12 -g transcript_id -a Homo_sapiens.GRCh38.96.chr.gtf -o countsinput.txt inputinfeccion_n1.bam inputinfeccion_n2.bam inputmock_n1.bam inputmock_n2.bam 2>> loginput.txt

+ b. Simplificar los counts
> cat countsinput.txt | cut -f 1,7-14 > simple_countsinput.txt
> cat countspulldown.txt | cut -f 1,7-14 > simple_countspulldown.txt

+ c. Llevar a cabo el análisis de expresión diferencial
> cat simple_countspulldown.txt | Rscript deseq1.r 2x2 > results_pulldown.txt 2>> log2.txt

+ d. Generar  Generate heatmap from the deseq1 normalized matrix
> cat norm-matrix-deseq1.txt | Rscript draw-heatmap.r > heatmap-norm-matrix-deseq1.pdf  2>> log2.txt

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/genes%20DE.png "Heatmap Input ")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/ModuloIII/genes%20deee.png "Heatmap PullDown")

***

