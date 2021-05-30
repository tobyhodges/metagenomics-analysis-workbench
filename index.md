---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---

{% include gh_variables.html %}

A lot of metagenomics analysis is done using command-line tools for three reasons:  
  
1) You will often be working with a large number of files, and working through the command-line rather than through a graphical user interface (GUI) allows you to automate repetitive tasks.

2) You will often need more compute power than is available on your personal computer, and connecting to and interacting with remote computers requires a command-line interface. 
  
3) You will often need to customize your analyses, and command-line tools often enable more customization than the corresponding GUI tools (if a GUI tool even exists).

In a previous lesson, you learned how to use the bash shell to interact with your computer through a command-line interface. In this lesson, you will be applying this new knowledge to 
carry out a common metagenomics workflow - identifying Operational Taxonomic Unities (OTUs) 
among samples taken from two metagenomes within a location. We will be starting with a set 
of sequenced reads (.fastq files), performing some quality control steps, assembly those 
reads into contigs, and ending by identifying and visualizing the OTUs among these samples.

As you progress through this lesson, keep in mind that, even if you aren’t going to be doing this same workflow in your research, you will be learning some very important lessons about using command-line bioinformatics tools. What you learn here will enable you to use a variety of bioinformatics tools with confidence and greatly enhance your research efficiency and productivity.

> ## Prerequisites
>
> This lesson assumes a working understanding of the bash shell. If you haven’t already 
> completed the [Shell Genomics lesson](https://nselem.github.io/shell-metagenomics/), and you aren’t 
> familiar with the bash shell, please review those materials before starting this lesson.
>
> This lesson also assumes some familiarity with biological concepts, 
> including the structure of DNA, nucleotide abbreviations, and the 
> concept microbiome.
> This lesson uses data hosted on an Amazon Machine Instance (AMI). Workshop participants
> will be given information on how to log-in to the AMI during the workshop. Learners using 
> these materials for self-directed study will need to set up their own AMI. Information 
> on setting up an AMI and accessing the required data is provided on the 
> [Metagenomics Workshop setup page](http://carpentries-incubator.github.io/metagenomics/setup.html).
>
{: .prereq}

> ## Things You Need To Know
>
> 0.  Stay calm, don't panic.
> 1.  Everything is going to be fine.
> 2.  We are learning together.
{: .checklist}

{% include links.md %}
