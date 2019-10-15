
# GITHUB
*Unidad 2 “Organización de un proyecto Bioinformático”*

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


##  Expresión Diferencial
*Unidad 7 “Análisis de Expresión diferencial a partir de secuencias de RNA”*

La Unidad 7 se basa en el análisis bioinformático de expresión diferencial de muestras de secuenciación de RNA con genoma de referencia, el cual permite determinar genes diferencialmente expresados sometidos a diferentes estímulos, permitiendo de esta manera interpretar de manera fácil mediante el análisis de gráficos y tablas. 

Para fines de este práctico, se utilziaron datos de 4 librerías de lecturas de arqueobacterias *Sulfolobus acidocaldarius*, el cual se le realizó un *knockdown* en el gen asociado a la formación de biopelículas. Las muestras fueron sometidas a diferentes ambientes, lo que posteriormente permitió diferenciar las distintas muestras:

* Wild type P
* Wild type B
* Mutant P
* Mutant B

### Pasos tutorial Unidad 7

#### 1. Crear variables
*Carpetas preexistentes que contienen la ubicación de las carpetas*

RAW=/shared/bioinfo1/common/raw_data/

ANN=/shared/bioinfo1/common/annot/

REF=/shared/bioinfo1/common/ref_genome/

#### 2. Crear carpetas de salida 
*Carpetas con rutas donde se almacenará información*

QC=../qc

FIL=../filtered

ALN=../alignment

CNT=../count

#### 3. Control de Calidad 
*Directorios donde se almacena información procesada. Luego se ejecuta el programa IlluQC_PRLL.pl, el cual genera un reporte completo de la calidad de las secuencias. La version PRLL permite ejecutar el comando usando distintos CPU al mismo tiempo. Una vez finalizado el análisis, el programa arroja distintos gráficos que representan la calidad de las lecturas, el contenido GC y otros datos necesarios para el análisis posterior de la secuenciación.* 

* Control de Calidad Muestra **Wild Type P**
> illuqc -se "$RAW/MW001_P.fastq" 5 A -onlystat -t 2 -o "wild_planctonic" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_QualRangePerBase.png "Cantidad de lecturas por base")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_avgQual.png "Valor promedio de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_gcDistribution.png "Distribucion de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_P.fastq_qualDistribution.png "Distribución de calidad")


***


* Control de Calidad Muestra **Wild Type B**
> illuqc -se "$RAW/MW001_B3.fastq" 5 A -onlystat -t 2 -o "wild_biofilm" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_QualRangePerBase.png "Cantidad de lecturas por base")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_gcDistribution.png "Distribución de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_qualDistribution.png "Distribución de calidad")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/MW001_B3.fastq_summary.png "Resumen de calidad de lecturas")


***


* Control de Calidad Muestra **Mutant P**
> illuqc -se "$RAW/0446_P.fastq" 5 A -onlystat -t 2 -o "mut_planctonic" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_QualRangePerBase.png "Cantidad de lecturas por base")

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_avgQual.png "Valor promedio de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_gcDistribution.png "Distribución de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_qualDistribution.png "Distribución de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_P.fastq_summary.png "Resumen de calidad de lecturas")


***


* Control de Calidad Muestra **Mutant B**
> illuqc -se "$RAW/0446_B3.fastq" 5 A -onlystat -t 2 -o "mut_biofilm" -c 10 &

![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_QualRangePerBase.png "Cantidad de lecturas por base")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_baseCompostion.png "Composición de Bases")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_gcDistribution.png "Distribución de Contenido de GC")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_qualDistribution.png "Distribución de calidad")


![alt text](https://github.com/mabayass/Tareas_Bioinfo2019_mby/blob/master/0446_B3.fastq_summary.png "Resumen de calidad de lecturas")


***
***

#### 4. Filtro de secuencias
*Luego de obtener los resultados del control de calidad de la secuenciación de RNA, las librerías son filtradas con el objetivo de eliminar lecturas con calidad menor de 20% en el 80% de la extensión, cuyos resultados genera librerías de lectura que seran utilizadas en el Alineamiento.*

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

***



#### 5. Alineamiento
*A partir de las librerias de lectura producidas en el punto anteior, se procede a hacer un alineamiento .*


