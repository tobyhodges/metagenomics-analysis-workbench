---
title: Reconstrucción de Genomas a partir de Metagenomas (MAGs)
teaching: 25
exercises: 15
questions: Cómo podemos obtener MAGs de buena calidad?
objectives: Tener una visión global sobre como reconstruir genomas de buena calidad
  a partir de metagenomas.
keypoints:
- Remover lecturas de hospedero mejora la calidad del metagenoma.
- La cobertura y composición nucleotídica son esneciales para el binning.
- Conocer bien nuestro objeto de estudio.
---

## Genomas a partir de metagenomas

La metagenómica hace referencia al estudio de todo el ADN de los organismos que se encuentran en un ambiente. La secuenciación de este material genético produce lecturas que pueden ensamblarse para conocer la diversidad microbiana y sus funciones.

Típicamente los metagenomas pueden estudiarse mediante dos aproximaciones:

- La clasificación taxonómica de contigs o lecturas y la inferencia metabólica de los contigs.
- La reconstrucción de genomas a a partir de metagenomas (MAGs), clasificación taxonómica y la inferencia metabólica de los MAGs.

En este apartado nos enfocaremos en la segunda aproximación. Los `MAGs` se reconstruyen a partir de un `ensamble metagenómico`,
los contigs de dicho ensamble se agrupan mediante la información de `cobertura y frecuencia de tetranucleótidos`.
Esta agrupación puede generar errores, por lo que es indispensable evaluar la calidad de los MAGs mediante la completitud
y redundancia de genes de copia única [MerenLab y col.](https://anvio.org/vocabulary/)

Para obtener MAGs podemos seguir el siguiente flujo de análisis:

<a href="fig/extrasMAGs/01.MAGs_workflow.png">
  <img src="fig/extrasMAGs/01.MAGs_workflow.png" alt="Flujo de trabajo para Metagenómica Centrada en Genomas" />
</a>

<a href="fig/extrasMAGs/01b.TNfCov.png">
  <img src="fig/extrasMAGs/01b.TNfCov.png" alt="Frecuencia de Tetranucleótidos y Profundidad" width="373" />
</a>
<br>

Ya que discutimos como seguir un flujo de análisis para reconstruir genomas entremos en acción, para ello analizaremos el metagenoma del pozol.

## El pozol

**El pozol** es un alimento ácido, fermentado a partir de maíz nixtamalizado, de importancia económica y cultural,
se consume desde tiempos prehispánicos y se ha estudiado desde los años 50s.

<a href="fig/extrasMAGs/02.Pozolhistoria.png">
  <img src="fig/extrasMAGs/02.Pozolhistoria.png" alt="Proceso de elaboración del pozol" />
</a>

<br>
Algunos puntos importantes que conocemos son:

<FONT COLOR="darkblue">
* No se inocula y al final de su fermentación tiene alta diversidad microbiana.<br>

- Es muy nutritivo, tiene un alto contenido de aminoácidos esenciales.<br>

- Es considerado como prebiótico, contiene fibras solubles y microorganismos benéficos para la salud intestinal humana.<br>

</FONT>

***

## Resolvamos preguntas biológicas mediante Metagenómica centrada en genomas

Imaginemos que se quiere impulsar la producción de esta bebida y para ello necesitan saber todo acerca de su naturaleza microbiana.

Una importante industria alimenticia los contacta como `expertos en ecología microbiana` y les pide ayuda para descubrir los siguientes puntos:
<br>

<FONT COLOR="darkblue">
* ¿Qué actores microbianos están presentes durante el proceso de fermentación?<br>
* ¿Cómo ocurre la bioconversión del maíz durante la fermentación, quién participa y cómo lo hace? ¿Qué funciones metabólicas están ocurriendo?<br>
* ¿Cambia la comunidad microbiana a lo largo del proceso?<br>
</FONT> <br>

La empresa secuenció cuatro puntos de fermentación de muestras que se obtuvieron en un mercado de Tabasco.
Las muestras se secuenciaron con Illumina NextSeq500 con lecturas pareadas de 75 pb.
Los datos están públicos bajo el Bioproject: [PRJNA648868](https://www.ebi.ac.uk/ena/browser/view/PRJNA648868)

<br>
<a href="fig/extrasMAGs/03.Pozol_fermentation.png">
  <img src="fig/extrasMAGs/03.Pozol_fermentation.png" alt="Puntos de fermentación" />
</a>

:::::::::::::::::::::::::::::::::::::::::  callout

## Limpieza de hospedero

Como las muestras contienen maíz, es indispensable remover las lecturas que correspondan a su genoma,
no hacerlo producirá un ensamble muy fragmentado, mayoritariamente del maíz y poco microbiano.
El autor del artículo amablemente nos proporcionó sus muestras libres del maíz y el código que usó
para ello está disponible en un repositorio público de [GitHub](https://github.com/RafaelLopez-Sanchez/pozol_shotgun).

[El artículo](https://www.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.001355): López-Sánchez et al., 2023. Analysing the dynamics of the bacterial community in pozol, a Mexican fermented corn dough.


::::::::::::::::::::::::::::::::::::::::::::::::::

***

## Espacio de trabajo

1. Entra a tu cuenta en el servidor y sitúate en tu `$HOME`
2. Obten los datos y la estructura de tu directorio del proyecto
3. Entra al directorio del proyecto

```bash
# ve al $HOME
cd

# descarga
#wget https://zenodo.org/records/13911654/files/taller_metagenomica_pozol.tar.gz?download=1 -O taller_metagenomica_pozol.tar.gz

# descomprime
#tar -xvzf taller_metagenomica_pozol.tar.gz

# Entra al directorio del proyecto
cd taller_metagenomica_pozol
```

Si en algún momento te pierdes entre directorios, puedes regresar al espacio principal asi:

:::::::::::::::::::::::::::::::::::::::::  callout

## Directorio principal del proyecto

```bash
cd && cd taller_metagenomica_pozol/
```

::::::::::::::::::::::::::::::::::::::::::::::::::

Cómo vamos a trabajar durante el taller?

:::::::::::::::::::::::::::::::::::::::::  callout

## Reglas del juego

- En este tutorial haremos el ejemplo corriendo la muestra de 48 hrs.
- Se formaran 6 equipos (2 de los tiempos 0, 9 y 24 hrs).
- Los equipos discutirán y presentarán sus resultados cuando se indique en el tutorial.
  

::::::::::::::::::::::::::::::::::::::::::::::::::

***

***

***

<p align="justify">

<FONT COLOR="darkblue">La presente práctica sólo es una representación del flujo de trabajo para el análisis metagenómico, sin embargo, `no sustituye los manuales` de cada programa y el flujo puede variar dependiendo del tipo de datos y pregunta de investigación.
De hecho para fines del taller, con frecuencia se utilizan las lineas de comando más simples para eficientar tiempo y recursos, tómalo en cuenta.</FONT>

</p>

Cada programa tiene una ayuda y un manual de usuario, es `importante` revisarlo y conocer cada parámetro que se ejecute. En terminal se puede consultar el manual con el comando `man` y también se puede consultar la ayuda con `-h` o `--help`, por ejemplo `fastqc -h`.

:::::::::::::::::::::::::::::::::::::::::  callout

## Para tenerlo presente

En bioinformática cualquier línea de comandos generará un resultado, de ahí a que esos resultados sean correctos puede haber una gran diferencia.
En cada paso detente a revisar la información de cada programa, lee el manual, visita foros de ayuda y selecciona los argumentos que se ajusten a las necesidades de tus datos.


::::::::::::::::::::::::::::::::::::::::::::::::::


