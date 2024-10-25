---
title: "Clasificación taxonómica y Metabolismo"
teaching: 1 h 
exercises: 2 h
questions:
- "Cuál es la identidad taxonómica de los MAGs?"
- "Qué funciones codifican?" 
objectives:
- "Obtener la identidad taxonómica de los MAGs y conocer su capacidad funcional"
keypoints:
- "La clasificación taxonómica por aproximación filogenética es altamente confiable"
- "Anotar los MAGs con diversas bases de datos mejora nuestro entendimiento de las funciones de los microorganismos"
---

## Taxonomía
### GTDB-tk

GTDB-Tk es una herramienta que asigna taxonomía a genomas utilizando la base de datos [GTDB](https://gtdb.ecogenomic.org/) (Genome Taxonomy Database).
Basado en árboles filogenéticos y medidas de ANI (Average Nucleotide Identity), GTDB-Tk clasifica genomas bacterianos y arqueanos, 
proporciona una taxonomía coherente y actualizada. Se utiliza mucho en el análisis de genomas y metagenomas.

<br>
<p style="text-align: center;">
  <a href="{{ page.root }}/fig/extrasMAGs/14.GTDB-tk_worflow.png">
    <img src="{{ page.root }}/fig/extrasMAGs/14.GTDB-tk_worflow.png" alt="Flujo de análisis de GTDB-tk"/>
  </a>
</p>
<br>

Recordemos que ya tenemos un set de *bins* refinados y desreplicados. Ahora vamos a asignarles identidad taxonómica, para ello vamos a correr GTDB-tk

<br>
Activa el ambiente
``` bash
conda activate gtdbtk-2.1.1
```
<br>
El directorio de resultados para gtdbtk ya lo tienes en tu carpeta de resultados. 
Para colocar los *bins* refinados y renombrados ejecuta el script `src/copiar_renombrarbins.sh`:

``` bash
 bash src/copiar_renombrarbins.sh
```
<br>
Ahora si, vamos a correr gtdbtk ...

``` bash
pip install numpy==1.19.5
gtdbtk classify_wf --genome_dir results/10.gtdbtk/bins/ --out_dir results/10.gtdbtk/ --cpus 4 -x fasta
```

No olvides desactivar el ambiente

``` bash
conda deactivate
```

> ## Resultado de GTDB-Tk
> Si gtdbtk está tomando mucho tiempo puedes parar el proceso con `ctrl + C` en tu teclado.
> El resultado final se encuentra en el directorio y archivo: `results/10.gtdbtk/gtdbtk.bac120.summary.tsv`.
{: .callout}

<br>
Después de ejecutar GTDB-tk, continuaremos en R para visualizar los datos.
Crea  un proyecto de R dentro del directorio `results/10.gtdbtk/`
<br>

``` r
library(tidyverse)
library(ggplot2)
# Leer la tabla ------------------------------------------------------------####
GTDBtk <- read.table("gtdbtk.bac120.summary.tsv", 
                    sep = "\t", header = TRUE, na.strings = "", 
                    stringsAsFactors = FALSE) %>%
  as_tibble()

# Transformar datos --------------------------------------------------------####

pozol_gtdbtk <- GTDBtk %>%
  select(user_genome, classification) %>%
  separate(classification, c(
    "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
    sep = ";") %>%
  rename(Bin_name = user_genome) %>%
  unite(Bin_name_2, c("Bin_name", "Phylum"), remove = FALSE) %>%
  select(Bin_name, Domain, Phylum, Class, Order, Family, Genus, Species)

# Guardamos los datos en un archivo de metadatos ---------------------------####
write.table(pozol_gtdbtk, file = "Metadatos.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Visualización de Datos ---------------------------------------------------####
GTDBtk_barplot <- pozol_gtdbtk %>%
  count(Phylum, Genus) %>%
  rename(Number_of_MAGs = n) %>%
  ggplot(aes(x = Phylum, y = Number_of_MAGs, fill = Genus)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal()

GTDBtk_barplot
```

> ## Discusión
> Qué otra estrategias implementarías para la clasificación taxonómica?
{: .callout}
<br>

> ## Ejercicio 3
> 
> Ahora te toca a tí.
> 
> * Reúnanse en equipos y obtengan la clasificación taxonómica de los MAGs que obtuvieron con la muestra que les tocó.
> * Discutan cada resultado obtenido.
> * En la [carpeta compartida de Drive](https://drive.google.com/drive/folders/1rg-zjuASg9D-goa2SlL3HXalqj3BQFNX) busquen la diapositiva para el Ejercicio 3.
> * En la diapositiva correspondiente resuman sus resultados obtenidos.
> 
> Tiempo de actividad (30 min)
> 
> Tiempo de presentación de resultados (2 min por equipo)
{: .challenge}

<br>
----------------------------------------------------------------------------------------------------------------

## Metabolismo

### PROKKA

[Prokka](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/annotation-with-prokka/slides-plain.html) es una herramienta útil, usa diferentes programas para predecir genes, secuencias codificantes, tRNAs, rRNAs. 
Hace la traducción de CDS a aminoácidos y asigna funciones usando varias bases de datos.


<br>
<p style="text-align: center;">
  <a href="{{ page.root }}/fig/extrasMAGs/15.Prokka_workflow.png">
    <img src="{{ page.root }}/fig/extrasMAGs/15.Prokka_workflow.png" alt="Prokka" width="673"/>
  </a>
</p>
<br>


Para correrlo vamos a activar el ambiente en el que se aloja.

``` bash
conda activate Prokka_Global
```
<br>
Tenemos el ambiente activo, ahora vamos a crear un directorio de resultados para prokka.

``` bash
mkdir -p results/11.prokka
```
<br>
Para correrlo, podemos hacer un ciclo que nos permita anotar todos los *bins.*

``` bash
nohup for FASTA in $(ls results/10.gtdbtk/bins/); do
    LOCUSTAG=$(basename $FASTA .fasta)
    prokka --locustag "${LOCUSTAG}_Scaffold"  --prefix $LOCUSTAG --addgenes --addmrn --cpus 4 --outdir "results/11.prokka/$LOCUSTAG" "results/10.gtdbtk/bins/$FASTA"
done &
```

> ## Explora
> 
> Mientras prokka se ejecuta en los bins que obtuviste, despliega la ayuda y discute:
> ¿ qué argumentos quitarías o agregarías?
>
> Cuáles te llamaron la atención?
{: .challenge}


<br>
Desactivemos el ambiente:

``` bash
conda deactivate
```

### KofamScan

Ahora que tenemos las proteínas predichas vamos a obtener más anotaciones útiles, usaremos kofam para esto.

[KofamScan](https://github.com/takaram/kofam_scan) es una herramienta de anotación, usa la base de datos KOfam de KEGG para obtener información sobre los genes que participan en diferentes rutas metabólicas.

Vamos a crear el directorio de resultados

``` bash
mkdir -p results/12.kofam
```


> ## Ejemplo de como correr KOfamScan
>
> KofamScan requiere mucho tiempo de ejecución.
> 
> Para efectos del taller nosotros **ya lo corrimos** y te proporcionaremos los resultados.
> 
> Pero te dejamos el bloque de código que usamos para este paso:
> 
> ``` bash
> for FAA in $(ls results/11.prokka/*/*.faa); do
>     name=$(basename $FAA .faa)
>     exec_annotation $FAA \
>         -o results/12.kofam/"$name.txt" \
>         --report-unannotated \
>         --cpu 4 \
>         --tmp-dir results/12.kofam/"tmp$name" \
>         -p /home/alumno1/kofam/db/profiles/ \
>         -k /home/alumno1/kofam/db/ko_list
> done
> # remover los directorios temporales
> #rm -r results/12.kofam/tmp*
> ```
{: .callout}

<br>
Estos resultados ya los tienes en el directorio `results/12.kofam`

<br>
<br>
Y ahora que ya tenemos los identificadores de KO para cada proteína, vamos a filtrar y graficar el metabolismo de los *bins*.

### RbiMs

[RbiMs](https://github.com/mirnavazquez/RbiMs) es un paquete de R muy útil para obtener la anotación de cada KEGG ID y generar plots de esta información. 
Puede trabajar con anotaciones de KOFAM, Interpro, PFAM o CAZY.


<br>
<p style="text-align: center;">
  <a href="{{ page.root }}/fig/extrasMAGs/16.rRbiMs-3.png">
    <img src="{{ page.root }}/fig/extrasMAGs/16.rRbiMs-3.png" alt="RbiMs"/>
  </a>
</p>
<br>

Vamos al editor de Rstudio para correr RbiMs ✨

``` r
library(tidyverse)
library(tidyr)
library(rbims)
library(readxl)

setwd("/home/alumnoX/taller_metagenomica_pozol/")

#A continuación, leemos los resultados de KEGG 
pozol_table <- read_ko(data_kofam = "results/12.kofam/") 

#y los mapeamos con la base de datos de KEGG:
pozol_mapp <- mapping_ko(pozol_table)

#Nos centraremos en las vías metabólicas relacionadas con la biosintesis de aminoacidos y vitaminas:

Overview <- c("Biosynthesis of amino acids", "Vitamin B6 metabolism")
Aminoacids_metabolism_pozol <- pozol_mapp %>%
  drop_na(Module_description) %>%
  get_subset_pathway(Pathway_description, Overview) 

#Visualizamos los datos con un gráfico de burbujas:

plot_bubble(tibble_ko = Aminoacids_metabolism_pozol,
            x_axis = Bin_name, 
            y_axis = Pathway_description,
            analysis = "KEGG",
            calc = "Percentage",
            range_size = c(1, 10),
            y_labs = FALSE,
            x_labs = FALSE)  

#Añadiremos metadatos, como la taxonomía:

Metadatos <- read_delim("results/10.gtdbtk/Metadatos.txt", delim = "\t")

#Y generaremos un gráfico de burbujas con metadatos:

plot_bubble(tibble_ko = Aminoacids_metabolism_pozol,
            x_axis = Bin_name, 
            y_axis = Pathway_description,
            analysis = "KEGG",
            data_experiment = Metadatos,
            calc = "Percentage",
            color_character = Family,
            range_size = c(1, 10),
            y_labs = FALSE,
            x_labs = FALSE) 

# Exploración de una Vía Específica
# podemos explorar una sola vía, como el “Secretion system,” y crear un mapa de calor para visualizar los genes relacionados con esta vía:

Secretion_system_pozol <- pozol_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(Cycle, "Secretion system")

#Y, finalmente, generamos un mapa de calor:

plot_heatmap(tibble_ko = Secretion_system_pozol, 
             y_axis = Genes,
             analysis = "KEGG",
             calc = "Binary")

#También podemos agregar metadatos para obtener una visión más completa:

plot_heatmap(tibble_ko = Secretion_system_pozol, 
             y_axis = Genes,
             data_experiment = Metadatos,
             order_x = Family,
             analysis = "KEGG",
             calc = "Binary")

plot_heatmap(tibble_ko = Secretion_system_pozol, 
             y_axis = Genes,
             data_experiment = Metadatos,
             order_y = Pathway_cycle,
             order_x = Family,
             analysis = "KEGG",
             calc = "Binary")

# Explorar
colnames(pozol_mapp) 

pozol_mapp %>%
  select(Cycle, Pathway_cycle, Pathway_description) %>%
  distinct()
```

<br>
<p style="text-align: center;">
  <a href="{{ page.root }}/fig/extrasMAGs/17.RBiMs_heatmap.png">
    <img src="{{ page.root }}/fig/extrasMAGs/17.RBiMs_heatmap.png" alt="RbiMs"/>
  </a>
</p>
<br>


### Antismash

Adicionalmente podrías anotar el metabolismo secundario de los *bins* siguiendo el flujo de análisis propuestos en la lección de [Minería Genómica](https://carpentries-incubator.github.io/genome-mining/02-antismash/index.html) de Software Carpentry.
