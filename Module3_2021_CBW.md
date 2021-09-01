---
layout: tutorial_page
permalink: /MIC_2021_Module3_lab
title: MIC 2021 Module 3 Lab
header1: Workshop Pages for Students
header2: Microbiome Analysis 2021
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/MIC_2021
description: MIC 2021 Module 3 Lab
---

This tutorial is part of the 2021 Canadian Bioinformatic Workshop.

**Authors**: Jacob Nearing and Morgan Langille

This tutorial is an updated version of the [2018 CBW Picrust2 tutorial](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2018-PICRUSt2-Tutorial/) original written by Gavin Douglas and Morgan Langille.

### Table of Contents
* [**Introduction**](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#introduction)
  * [Teaching objectives](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#teaching-objectives)
  * [Bioinformatic tool citations](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#bioinformatic-tool-citations)
  * [Activating the conda environment](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#activating-the-conda-environment)
  * [Exploring input files](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#exploring-input-files)
* [**Read placement**](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#read-placement)
* [**Hidden-state prediction**](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#hidden-state-prediction)
* [**Metagenome prediction**](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#metagenome-pipeline)
* [**Pathway-level inference**](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#pathway-inference)
* [**Add Descriptions**](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial#add-descriptions)


## Introduction
PICRUSt2 is a tool that predicts the abundance of gene families and higher-level pathways present in a microbial community based on amplicon sequencing data. This tutorial will demonstrate each major step in the PICRUSt2 pipeline. Throughout the tutorial there will be questions that are meant to aid understanding. [The answers are on this page](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial-Answers), which you might want to open in a different tab.

For this lab we will be processing 16S sequencing of ileum and rectum biopsy samples from the [Inflammatory Bowel Disease Multi'omics Database](https://ibdmdb.org/). Although shotgun metagenomics is becoming cheaper it typically isn't done with biopsy samples since it requires an enormous amount of sequencing (due to the large amount of host DNA!). 16S sequencing doesn't suffer from this issue since only the 16S is amplified - there typically is no host contamination.

This tutorial was meant for PICRUSt v2.4.1.
### Teaching objectives

1. Run and understand major steps of PICRUSt2 workflow.
2. Understand what files PICRUSt2 outputs and what the abundances in each one refer to.

### Bioinformatic tool citations
* **castor** ([website](https://cran.r-project.org/web/packages/castor/index.html), [paper](https://academic.oup.com/bioinformatics/article-abstract/34/6/1053/4582279?redirectedFrom=fulltext)) 
* **EPA-NG** ([website](https://github.com/Pbdas/epa-ng), [pre-print](https://www.biorxiv.org/content/early/2018/03/29/291658))
* **gappa** ([website](https://github.com/lczech/gappa)) - wrapper for the [genesis library](https://github.com/lczech/genesis).
* **MinPath** ([website](http://omics.informatics.indiana.edu/MinPath/), [paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000465))
* **PaPaRa** ([website](https://github.com/sim82/papara_nt), [paper](https://academic.oup.com/bioinformatics/article/27/15/2068/400617))
* **PICRUSt2** ([website](https://github.com/picrust/picrust2/wiki), [PICRUSt2 paper](https://www.nature.com/articles/s41587-020-0548-6))
* **HMMER** ([website](http://hmmer.org/), [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195))

### Activating the conda environment

PICRUSt2 is already installed on each of the students instances in a [conda](https://conda.io/docs/user-guide/getting-started.html) environment, which is a tool for managing environments. This is helpful since sometimes different tools require different versions of packages as dependencies. It also allows users to easily install software with a simple command. The environment we'll be using for this tutorial is called ```picrust2```.

You can activate the environment for this tutorial with this command:
```
conda activate picrust2
```

When you're **finished** this tutorial you can deactivate the environment with this command:

```
conda deactivate
```


### Exploring input files

Two input files are required to run PICRUSt2: (1) a FASTA file of each study sequence and (2) the abundance of each study sequence in each sample. The workflow will work with _any_ sequences that can be reliably placed in the reference tree (which is based on 16S rRNA genes [16S] by default), but typically users will input either the representative sequences of OTUs or amplicon sequence variants (ASVs). Think back to your earlier sessions on 16S rRNA sequencing to remind yourself what ASVs/OTUs are!


The first step will be to make a new directory for this tutorial so that we can keep all of our files organized. We will also download the required files for this tutorial using the wget command.
```
mkdir ~/workspace/picrust2_tut
cd ~/workspace/picrust2_tut
mkdir input_files
wget https://www.dropbox.com/s/cjqy61z4g1pylbz/ASV_abun.tsv?dl=1 -O input_files/ASV_abun.tsv
wget https://www.dropbox.com/s/j865g782dc3t96g/ASVs.fna?dl=1 -O input_files/ASVs.fna
wget https://www.dropbox.com/s/mhly9pf53vsfx7l/picrust2_lab_metadata.tsv?dl=1 -O input_files/picrust2_lab_metadata.tsv
```
The above wget commands will have downloaded three files.

The first file ```ASV_abun.tsv``` represents each sample as a column (sample ids are in first row), each ASV as a row (ASV ids are in first column), and the numbers in the table as counts of ASVs per sample.

View the ```ASV_abun.tsv``` file with the ```less``` command and the option ```-S``` which turns off line wrapping (use the left and right arrow keys to scroll to see the entire file):

```
less -S input_files/ASV_abun.tsv
```

You can also see how many lines there in the sequence table with this command:

```
wc -l input_files/ASV_abun.tsv
```

There are 100 sequences in this table (the first line is the header!)


The second file is named ```picrust2_lab_metadata.tsv```. This file contains all the metadata for the samples of interest. We will use this file along with some pre-made Rscripts to visualize how our predictions differ between Rectum and Ileum samples from individuals with CD and nonIBD.

The final file is named ```ASVs.fna```. This file contains all of the representative sequences for the ASVs in our ASV abundance table. You will notice that this file ends in ```.fna``` and not ```.fasta```. You may encounter this on your own but do not worry as ```.fna``` format is just another file identifier for FASTA formatted files.

**Question 1: How many samples are in this dataset?**

**Question 2: How could you check how many sequences are in the FASTA file ```input_files/ASVs.fna```?**

## Read placement

The first step of the workflow involves taking the study sequences and putting them into the pre-existing phylogenetic tree based on the 16S sequences from the reference genomes with the script ```place_seqs.py```. As with all of the scripts in this workflow you can see the arguments and options for this tool by typing ```place_seqs.py --help```. The default reference tree is based on full-length 16S sequences parsed from genomes in the [Integrated Microbial Genomes database](https://img.jgi.doe.gov/). The multiple-sequence alignment (MSA) of these reference 16S sequences is also required as input since it is required for the first placement step (done by _HMMER_). However, luckily for us this script automatically points to this MSA file that was included with our PICRUSt2 installation!


To run the sequence placement algorithm and produce a new phylogenetic tree with our sequences placed within them we can run the following command.

```
place_seqs.py -s input_files/ASVs.fna -o tutorial.tre -p 4 --verbose
```

These are the inputs to this tool (besides changing the default reference files we talked about eariler):
* ```-s FASTA```: your study sequences (i.e. FASTA of amplicon sequence variants or operational taxonomic units)
* ```-o TREEFILE```: Output tree with placed study sequences.
* ```-p INT```: Number of threads to use.
* ```--verbose```: Option to specify that wrapped commands will be printed to screen (useful for troubleshooting!).

The output of this tool is the tree ```tutorial.tre```. Take a look at this file with ```less```. The tree is in [newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html), which captures how individual tips on the tree are related topologically and also the branch-length distances between them. Importantly, ```place_seqs.py``` takes the **most likely** placement for each input sequence, but that doesn't mean it's correct!

## Hidden-state prediction

Now that we have the study sequences placed in the reference tree we can proceed to the heart of the PICRUSt2 pipeline: hidden-state prediction (HSP). For each placed study sequence the predicted abundances of gene families of interest will be predicted, which below will be Enzyme Classification (E.C.) numbers.

We will also need to predict how many 16S copies are expected to be in the genome corresponding to each study sequence. Predicting how many 16S copies there are is important since this will be used in the next step to normalize the abundance of each sequence. We will also calculate the nearest-sequenced taxon index (NSTI), which is the measure of how distant each study sequence from the nearest reference sequence in the tree.

This command will run HSP using maximum parsimony to predict the number of 16S genes per predicted genome and will calculate NSTI values for each study sequence:

```
hsp.py -i 16S -t tutorial.tre -o 16S_predicted.tsv -m mp -p 4 -n
```

Similarly, this command will generate predictions on how many copies of each gene family are found in each predicted genome. We wont bother calculating NSTI values again since it would be the same as above.

```
hsp.py -i EC -t tutorial.tre -o EC_predicted.tsv -p 4 -m mp
```

The input arguments and options are:
* ```-t TREEFILE```: Newick tree with study sequences placed amongst reference sequences.
* ```-i TRAIT_OPTION```: Which default pre-calculated count table to use (one of '16S', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM', 'PHENO')
* ```-o FILENAME```: Named of output file containing predicted counts.
* ```-m METHOD```: Hidden-state prediction method to use.
* ```-p INT```: Number of processes to run in parallel.
* ```-n```: Indicates that Nearest-sequenced taxon index (NSTI) values should be calculated.

There are a few other options which you can look at using the command ```hsp.py --help``` however I would suggest users stick to the default values for the other arguments unless they are using a custom database.

Next, it's a good idea to check whether any of your input sequences have very high NSTI values. It's possible that sequences with high NSTI values (e.g. > 1) could be uncharacterized taxa, but often they are simply garbage sequences.

Lets first inspect the ```16S_predicted.tsv``` file using the ```less``` command.

```
less 16S_predicted.tsv
```

By doing this you will notice that there are three columns within the file: sequence, 16S_rRNA_Count and metadata_NSTI. You could sort this table by the NSTI values in the 3rd column using the sort command:

```
sort -k 3 16S_predicted.tsv
```

**Question 3: What sequence has the highest NSTI value?**

**Question 4: Which sequence is predicted to have the most copies of the 16S gene? How many copies does it have?**

You can also look at the predictions for ECs:

```
less -S EC_predicted.tsv
```

Each EC number is a different column and thus there are many more columns compared to the 16S copy number predictions!

## Metagenome pipeline
Now that we have the predicted _genomes_ for each study sequence we can predict the _metagenomes_ for each sample. The basic approach here is to multiply the predicted number of genes by the abundance of each corresponding "genome" (i.e. the abundance of the study sequence). However, before predicting the metagenomes the sequence abundances are normalized by the number of marker gene copies in each predicted genome (i.e. how many 16S copies they have).

```
metagenome_pipeline.py -i input_files/ASV_abun.tsv -m 16S_predicted.tsv -f EC_predicted.tsv -o EC_metagenome_out --strat_out
```

The input arguments/options are:
* ```-i STUDY.tsv``` - Table of ASV abundances across all samples (either TSV or BIOM format).
* ```-m MARKER_PREDICTED.txt``` - Output predicted 16S copy numbers (or other marker) for all study sequences.
* ```-f FUNC_PREDICTED.txt``` - Output predicted functional abundances for all study sequences.
* ```--max_nsti INT``` - Max NSTI values per study sequence (sequences with values above this cut-off will be removed). Default: 2.
* ```-o OUTPUT_DIR``` - Output directory.
* ```--strat_out``` - Indicates we would like to generate an additional file with the functional contributions from each individual ASV.

The above command output the predicted metagenomes to the folder ```EC_metagenome_out``` by default. There are 4 output files in this directory:
* ```pred_metagenome_strat.tsv.gz```: Predicted functional contribution to the metagenome from each individual ASV
* ```pred_metagenome_unstrat.tsv.gz```: predicted metagenomes over all sequences.
* ```seqtab_norm.tsv.gz```: sequence table normalized by predicted number of marker genes.
* ```weighted_nsti.tsv.gz```: weighted NSTI values per sample (over all sequences weighted by the abundance in that sample).

You will notice that these files are compressed in gzip format. We will go ahead and uncompress them with the following command:
```
gunzip EC_metagenome_out/*
```

Take a look at the normalized sequence abundance table with ```less - S EC_metagenome_out/seqtab_norm.tsv```. The abundances in this file are the number of reads divided by the predicted number of marker genes.

**Question 5: What is the abundance of sequence ```132b62bcbd3f280208e314eeab566a4d``` in sample ```HMP.2.206659``` in this normalized table? Given that the unnormalized abundance was 11 how many predicted 16S sequences must there have been?**

The output weighted NSTI values are also useful to explore.

You will also notice that this part of the pipeline output two different predicted metagenome files. The first file ```pred_metagenome_unstrat.tsv``` represents a table similar to what we saw with the ASV_abun.tsv file. On columns we have each different sample and on the rows we have the predicted normalized abundance for each enzyme classification ID. This is one of the major output files you would use to begin looking for differences in functions across your samples. 

```
less -S EC_metagenome_out/pred_metagenome_unstrat.tsv
```

**Question 6: Why would the column sums (total predicted abundance in each sample) in the unstratified predictions not be equal even if the input sequence abundance table was rarified so that all samples had the same depth?**

The second predicted metagenome file ```pred_metagenome_strat.tsv``` outputs the breakdown of which taxa are contributing to each predicted function (e.g. EC number). This is useful to determine which taxa are causing the shift in a particular function.

## Pathway inference

The final step in the PICRUSt2 pipeline is to infer pathway abundances based on the presence of gene families. This step is done by wrapping [MinPath](http://omics.informatics.indiana.edu/MinPath/), which determines the minimum number of pathways that can be explained given the presence of gene families. The default is to map the EC numbers to Metacyc reactions and then to Metacyc Pathways.

You can run this step on the Stratified table using this command:

```
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv -o pathways_out -p 4
```

The inputs to this tool are:
* ```-i EC_metagenome_out/pred_metagenome_contrib.tsv``` - Stratified output of ```metagenome_pipeline.py``` step. (Note you can also input the unstratified table, however if a stratified table is input then it will generate both stratified and unstratified pathway tables)
* ```-o pathways_out``` - The name of the folder to output files into.
* ```-p INT``` - Number of processes to run simultaneously.


Once again both stratified and unstratified output files are created: ```pathways_out/path_abun_contrib.tsv.gz``` and ```pathways_out/path_abun_unstrat.tsv.gz``` respectively.

The abundance of pathways in these files is based on the mean abundance of the required gene families.

You can decompress these files and view them using previously used commands (e.g. ```gunzip``` and ```less```)

## Add Descriptions

The default outputs from PICRUSt are simple functional ids which often do not provide much information by themselves. You can manually look up these ids on their website (e.g. KEGG, Metacyc, etc.). However, you can also add a brief description to these files using the `add_descriptions.py` script which will add a column to your table of gene family or pathway abundances corresponding to a quick description of each functional category. These descriptions are in picrust2/default_files/description_mapfiles. You can also use custom mapping files.

To add a description column to E.C. number and MetaCyc pathway abundance tables respectively:

```
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv
```

```
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv
```

Arguments/options:

* ```-i: Input function abundance table``` - Can be stratified or unstratified. Can be compressed or uncompressed.
* ```-o: Output function abundance table with added description column.``` - Will be compressed if .gz is given
* ```-m: Default mapping table to use``` - Options are METACYC,COG,EC,KO,PFAM,TIGRFAM

**Question 7: What is the description for the pathway id PWY-6317?**

These tables can then be loaded into various visualization and statistical programs for further analysis. Basic browsing and visualization could be done in a spreadsheet program like Excel, while more advanced statistical and visualization could be conducted in R. 



