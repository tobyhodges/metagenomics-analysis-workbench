---
title: "Taxonomic Assignation"
teaching: 30
exercises: 15
questions:
- "How can I assign a taxonomy to my contigs?"
objectives:
- "Understand how taxonomic assignation works."
keypoints:
- "A database with previous gathered knowledge (genomes) is needed for taxonomic assignation."
- "Kraken is a program for taxonomic assignation."
---
## What is taxonomic assignation?

Taxonomic assignation is the process of assigning an Operational Taxonomic
Unit (OTUs, that is, groups of related individuals) to sequences, that can be 
reads or contigs. To assign an OTU to a sequence it is compared against a database, 
but this comparison can be done in different ways. The comparison database in 
this assignation process must be constructed using 
complete genomes. There are many programs for doing taxonomic mapping, 
almost all of them follows one of the next strategies:  

1. BLAST: Using BLAST or DIAMOND, these mappers search for the most likely hit 
for each sequence within a database of genomes (i.e. mapping). This strategy is slow.    
  
2. K-mers: A genome database is broken into pieces of length k, so as to be able to 
search for unique pieces by taxonomic group, from lowest common ancestor (LCA), 
passing through phylum to species. Then, the algorithm 
breakes the query sequence (reads, contigs) into pieces of length k,
look for where these are placed within the tree and make the 
classification with the most probable position.  

<a href="{{ page.root }}/fig/03-07-01.png">
  <img src="{{ page.root }}/fig/03-07-01" alt="Cog Metagenome" />
</a>
<!-- ![image](https://user-images.githubusercontent.com/67386612/119909687-20c90800-bf1b-11eb-8878-23b3773016ff.png)!-->
###### Figure 1. Lowest common ancestor assignation example.

3. Markers: They look for markers of a database made a priori in the sequences 
to be classified and assign the taxonomy depending on the hits obtained.    

A key result when you do taxonomic assignation of metagenomes is the abundance of each taxa or OTU in your sample. The absolute abundance of a taxon is the number of sequences (reads or contigs, depending on what you did) assigned to it. And its relative abundance is the proportion of sequences assigned to it. It is important to be aware of the many biases that that can skew the abundances along the metagenomics workflow, shown in the figure, and that because of them we may not be obtaining the real abundance of the organisms in the sample.

<a href="{{ page.root }}/fig/03-07-02.png">
  <img src="{{ page.root }}/fig/03-07-02.png" alt="Cog Metagenome" />
</a>
###### Figure 2. Abundance biases during a metagenomics protocol

## Using Kraken 2

[Kraken 2](https://ccb.jhu.edu/software/kraken2/) is the newest version of Kraken, 
a taxonomic classification system using exact k-mer matches to achieve 
high accuracy and fast classification speeds. `kraken2` is already installed in the metagenome
environment, lets have a look at `kraken2` help.  
 
~~~  
$ kraken2  
~~~ 
{: .bash}

~~~
Need to specify input filenames!
Usage: kraken2 [options] <filename(s)>

Options:
  --db NAME               Name for Kraken 2 DB
                          (default: none)
  --threads NUM           Number of threads (default: 1)
  --quick                 Quick operation (use first hit or hits)
  --unclassified-out FILENAME
                          Print unclassified sequences to filename
  --classified-out FILENAME
                          Print classified sequences to filename
  --output FILENAME       Print output to filename (default: stdout); "-" will
                          suppress normal output
  --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                          in [0, 1].
  --minimum-base-quality NUM
                          Minimum base quality used in classification (def: 0,
                          only effective with FASTQ input).
  --report FILENAME       Print a report with aggregrate counts/clade to file
  --use-mpa-style         With --report, format report output like Kraken 1's
                          kraken-mpa-report
  --report-zero-counts    With --report, report counts for ALL taxa, even if
                          counts are zero
  --report-minimizer-data With --report, report minimizer and distinct minimizer
                          count information in addition to normal Kraken report
  --memory-mapping        Avoids loading database into RAM
  --paired                The filenames provided have paired-end reads
  --use-names             Print scientific names instead of just taxids
  --gzip-compressed       Input files are compressed with gzip
  --bzip2-compressed      Input files are compressed with bzip2
  --minimum-hit-groups NUM
                          Minimum number of hit groups (overlapping k-mers
                          sharing the same minimizer) needed to make a call
                          (default: 2)
  --help                  Print this message
  --version               Print version information

If none of the *-compressed flags are specified, and the filename provided
is a regular file, automatic format detection is attempted.

~~~  
{: .output}

In addition to our input files we also need a database with which to compare them. There are [several databases](http://ccb.jhu.edu/software/kraken2/downloads.shtml) 
compatible to be used with kraken2 in the taxonomical assignation process. 

> ## Very important to know your database!
> The database you use will determine the result you get for your data.
> Imagine you are searching for a lineage that was recently discovered and it is not part of the available databases. Would you find it?
{: .callout}

Minikraken is a popular database that attempts to conserve its sensitivity 
despite its small size (Needs 8GB of RAM for the assignation). Unfortunately although it is much smaller that most databases, it is not small enough to be run by the machines we are using, so we won't be able to run `kraken2`. We can check our available RAM with `free -h`to be shure of this.
~~~
$ free -h
~~~
{: .bash}

~~~
              total        used        free      shared  buff/cache   available
Mem:           3.9G        272M        3.3G         48M        251M        3.3G
Swap:            0B          0B          0B
~~~
{: .output}

If we were to download the database we would use the following command:
~~~
$ curl -O ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz         
$ tar -xvzf minikraken2_v2_8GB_201904.tgz 
~~~
{: .do not run this}

> ## Exercise 1: Remembering commands
> 
> What is the command `tar` doing to the file `minikraken2_v2_8GB_201904.tgz`.  
> 
>> ## Solution
>> `tar` command is used in linux to decompress files, so in this case it would
>> extract the content of the compressed file  `minikraken2_v2_8GB_201904.tgz`  
>> 
> {: .solution}
{: .challenge}                             
                             
## Taxonomic assignation of metagenomic reads

As we have learned, taxonomic assignation can be attempted before the assembly process. 
In this case we would use FASTQ files as inputs, which would be 
`JP4D_R1.trim.fastq.gz` and `JP4D_R2.trim.fastq.gz`. And the outputs would be two files: the report
`JP4D.report` and the kraken file `JP4D.kraken`.  
  
To run kraken2 we would use a command like this:  
~~~
$ mkdir TAXONOMY_READS
$ kraken2 --db kraken-db --threads 8 --paired --fastq-input JP4D_R1.trim.fastq.gz JP4D_R2.trim.fastq.gz --o
utput TAXONOMY_READS/JP4D.kraken --report TAXONOMY_READS/JP4D.report
~~~
{: .do not run this}

Since we can't run `kraken2` here, we precomputed its results in a server, i.e. a more powerful machine. 
In the server, after we assembled the metagenome for this sample, we ran `kraken2` and obtained`JP4D-krakne.kraken` and `JP4D.report`.

Let's look at the precomputed outputs of `kraken2` in our assembled metagenome.  
~~~
head ~/dc_workshop/taxonomy/JP4D.kraken  
~~~
{: .bash}

~~~
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:19691:2037	0	250|251	0:216 |:| 0:217
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:14127:2052	0	250|238	0:216 |:| 0:204
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:14766:2063	0	251|251	0:217 |:| 0:217
C	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:15697:2078	2219696	250|120	0:28 350054:5 1224:2 0:1 2:5 0:77 2219696:5 0:93 |:| 379:4 0:82
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:15529:2080	0	250|149	0:216 |:| 0:115
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:14172:2086	0	251|250	0:217 |:| 0:216
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:17552:2088	0	251|249	0:217 |:| 0:215
U	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:14217:2104	0	251|227	0:217 |:| 0:193
C	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:15110:2108	2109625	136|169	0:51 31989:5 2109625:7 0:39 |:| 0:5 74033:2 31989:5 1077935:1 31989:7 0:7 60890:2 0:105 2109625:1
C	MISEQ-LAB244-W7:156:000000000-A80CV:1:1101:19558:2111	119045	251|133	0:18 1224:9 2:5 119045:4 0:181 |:| 0:99
~~~
{: .output}

As we can see, the kraken file is not very readable. So let's look at the report file:

> ## Reading a Kraken report
>
> 1. Percentage of reads covered by the clade rooted at this taxon
> 2. Number of reads covered by the clade rooted at this taxon
> 3. Number of reads assigned directly to this taxon
> 4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
> 5. NCBI taxonomy ID
> 6. Indented scientific name
{: .callout}

~~~
head ~/dc_workshop/taxonomy/JP4D.report
~~~
{: .bash} 
~~~
 78.13	587119	587119	U	0	unclassified
 21.87	164308	1166	R	1	root
 21.64	162584	0	R1	131567	  cellular organisms
 21.64	162584	3225	D	2	    Bacteria
 18.21	136871	3411	P	1224	      Proteobacteria
 14.21	106746	3663	C	28211	        Alphaproteobacteria
  7.71	57950	21	O	204455	          Rhodobacterales
  7.66	57527	6551	F	31989	            Rhodobacteraceae
  1.23	9235	420	G	1060	              Rhodobacter
  0.76	5733	4446	S	1063	                Rhodobacter sphaeroides
~~~
{: .output}  

## Taxonomic assignation of the contigs of a MAG
 
We now have the taxonomic identity of the reads of the whole metagenome, but we do not know to which taxon our MAGs correspond to. For this we have to make the taxonomic assignation with their contigs instead of its reads, because we do not have the reads that correspond to a MAG separated from the reads of the entire sample. 

For this the `kraken2` is a little bit different; here we can look at the command for the `JP4D.001.fasta` MAG:
  
~~~
$ mkdir TAXONOMY_MAG
$ kraken2 --db kraken-db --threads 12 -input JP4D.001.fasta --output TAXONOMY_MAG/JP4D.001.kraken --report TAXONOMY_MAG/JP4D.001.report
~~~
{: .do not run this}

The results of this are pre-computed in the `~/dc_workshop/taxonomy/mags_taxonomy/` directory
~~~
$ cd ~/dc_workshop/taxonomy/mags_taxonomy
$ ls
~~~
{: .bash} 
~~~
JP4D.001.kraken
JP4D.001.report
~~~
{: .output} 

~~~
more ~/dc_workshop/taxonomy/mags_taxonomy/JP4D.001.report
~~~
{: .bash} 
~~~
 50.96	955	955	U	0	unclassified
 49.04	919	1	R	1	root
 48.83	915	0	R1	131567	  cellular organisms
 48.83	915	16	D	2	    Bacteria
 44.40	832	52	P	1224	      Proteobacteria
 19.37	363	16	C	28216	        Betaproteobacteria
 16.22	304	17	O	80840	          Burkholderiales
  5.66	106	12	F	506	            Alcaligenaceae
  2.72	51	3	G	517	              Bordetella
  1.12	21	21	S	2163011	                Bordetella sp. HZ20
  .
  .
  .
~~~
{: .output}  

By looking at the report, we can see that half of the contigs are unclassified, and that a very little proportion of contigs have been assigned an OTU. This is weird because we expected to have only one genome in the bin.

Just to exemplify how a report of a complete and not contaminated MAG should look like, let's look at the report of this MAG from another study:
~~~
100.00	108	0	R	1	root
100.00	108	0	R1	131567	  cellular organisms
100.00	108	0	D	2	    Bacteria
100.00	108	0	P	1224	      Proteobacteria
100.00	108	0	C	28211	        Alphaproteobacteria
100.00	108	0	O	356	          Rhizobiales
100.00	108	0	F	41294	            Bradyrhizobiaceae
100.00	108	0	G	374	              Bradyrhizobium
100.00	108	108	S	2057741	                Bradyrhizobium sp. SK17
~~~
{: .output} 

> ## Discussion: 
>
> Why do you think we found so many OTUs in this bin? 
{: .discussion}

## Visualization of taxonomic assignation results  
We have reached the tsv files, the final step in our metagenomic pipeline showed in [lesson-3](https://carpentries-incubator.github.io/metagenomics/03-assessing-read-quality/index.html).  
After we have the taxonomy assignation what follows is some visualization of our results. 
[Krona](https://github.com/marbl/Krona/wiki) is a hierarchical data visualization software. Krona allows data to be explored with zooming, multi-layered pie charts and includes support for several bioinformatics tools and raw data formats. To use Krona in our results, lets go first into our taxonomy directory, which contains the precalculated Kraken outputs.  

### Krona  
~~~
$ cd ~/dc_workshop/taxonomy/mags_taxonomy
~~~
{: .bash}  

Krona is called with the `ktImportTaxonomy` command that needs an input and an output file.  
In our case we will create the input file with the columns three and four from `JP4D.001.kraken` file.     
~~~
$ cut -f2,3 JP4D.001.kraken > JP4D.001.krona.input
~~~
{: .language-bash}  

Now we call Krona in our `JP4D.001.krona.input` file and save results in `JP4D.001.krona.out.html`.  
~~~
$ ktImportTaxonomy JP4D.001.krona.input -o JP4D.001.krona.out.html
~~~
{: .language-bash}  

Instead of the output we expected and error appeared. 
~~~
Loading taxonomy...
Taxonomy not found in /home/dcuser/.miniconda3/envs/metagenomics/opt/krona/taxonomy. Was updateTaxonomy.sh run? at /home/dcuser/.miniconda3/envs/metagenomics/opt/krona/scripts/../lib/KronaTools.pm line 1540.
~~~
{: .error}  

It seems that a necessary command for Krona to work was not executed, so let's do that. But we need to deactivate our environment first.

~~~
$ conda deactivate
$ bash /home/dcuser/.miniconda3/envs/metagenomics/opt/krona/updateTaxonomy.sh 
~~~
{: .language-bash}  

Once it's done we activate the environment and try again.
~~~
$ conda activate metagenomics
$ ktImportTaxonomy JP4D.001.krona.input -o JP4D.001.krona.out.html
~~~
{: .language-bash}  

~~~
Loading taxonomy...
Importing JP4D.001.krona.input...
   [ WARNING ]  The following taxonomy IDs were not found in the local database and were set to root
                (if they were recently added to NCBI, use updateTaxonomy.sh to update the local
                database): 1804984 2109625 2259134
~~~
{: .output}  

And finally, open another terminal in your local computer,download the Krona output and open it on a browser and explore this visualization tool.
~~~
$ scp dcuser@ec2-3-235-238-92.compute-1.amazonaws.com:~/dc_workshop/taxonomy/JP4D.001.krona.out.html . 
~~~
{: .bash}  

<a href="{{ page.root }}/fig/03-07-03.png">
  <img src="{{ page.root }}/fig/03-07-03.png" alt="Krona Visualization" />
</a>

### Pavian
Pavian is another visualization tool that allows comparison between multiple samples. 
Pavian should be locally installed and needs R and Shiny, 
but we can try the [Pavian demo WebSite](https://fbreitwieser.shinyapps.io/pavian/) 
to visualize our results.  

First we need to download the files needed as inputs in Pavian, this time we will visualize the assignation of the reads of both samples:
`JC1A.report` and `JP4D.report`.  
This files corresponds to our Kraken reports. Again in our local machine lets use `scp` command.  
~~~
$ scp dcuser@ec2-3-235-238-92.compute-1.amazonaws.com:~/dc_workshop/report/*report . 
~~~
{: .language-bash}

We go to the [Pavian demo WebSite](https://fbreitwieser.shinyapps.io/pavian/), click on Browse and choose our reports.

<a href="{{ page.root }}/fig/03-07-04.PNG">
  <img src="{{ page.root }}/fig/03-07-04.PNG" alt="upload Pavian" />
</a>

We click on the Results Overview tab.

<a href="{{ page.root }}/fig/03-07-05.PNG">
  <img src="{{ page.root }}/fig/03-07-05.PNG" alt="Results Overview" />
</a>

We click on the Sample tab.

<a href="{{ page.root }}/fig/03-07-06.PNG">
  <img src="{{ page.root }}/fig/03-07-06.PNG" alt="sample" />
</a>

We can look at the abundance of a specific taxon by clicking on it.

<a href="{{ page.root }}/fig/03-07-07.PNG">
  <img src="{{ page.root }}/fig/03-07-07.PNG" alt="Sample Selected" />
</a>

We can look at a comparison of both our samples in the Comparison tab. 

<a href="{{ page.root }}/fig/03-07-08.PNG">
  <img src="{{ page.root }}/fig/03-07-08.PNG" alt="Comparison" />
</a>


> ## Discussion: Taxonomic level of assignation
>
> What do you think is harder to assign, a species (like _E. coli_) or a phylum (like Proteobacteria)?
{: .discussion}

                             
{% include links.md %}
