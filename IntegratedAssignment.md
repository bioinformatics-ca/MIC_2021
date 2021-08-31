# CBW Microbiome Analysis 2021

# Integrated Assignment

## **Developed by:** Justin Jia & Muhammad Zohaib Anwar

# Introduction

By now you’ve been introduced to a few pipelines for analyzing both 16S amplicon and shotgun metagenomics sequencing data. For this assignment we’ll be re-visiting these pipelines, but we will not be providing all of the commands for you to copy-and-paste. Two key skills for bioinformaticians to have is to be able to adapt existing pipelines on your data and to quickly learn how to use new tools. You’ll need to showcase both of these skills to complete this assignment!

For your integrated assignment, we will use our knowledge gained from this workshop to explore a (sub-sampled) dataset collected from ***Shaiber. A et al.***

[Functional and genetic markers of niche partitioning among enigmatic members of the human oral microbiome - Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w)

## A quick introduction to this dataset

Microbiologists often studies the ecology of oral microbes in the Human Microbiome Project (HMP). They found extensive site specificity among oral microbes, making it a fascinating environment to study microbial colonization. The oral cavity is a rich environment with multiple distinct niches in a relatively small space partially due to its diverse anatomy, hard and soft tissue structures, and exposure to exogenous factors. Oral microbes form complex communities that show remarkable patterns of horizontal and vertical transmission across humans and animals. For details, please review the manuscript above.

*Shaiber. A et al.*  sequenced amplicons (**16S marker gene**), **shotgun metagenomics data and long read genomes.** We will be using 8 samples from 16S rRNA dataset and shotgun metagenomes in this assignment. Table below gives an insight into selected samples.

[Samples and files for Assignment ](https://www.notion.so/327e5e3fb1d149419655d1ed08936f4e)

For shotgun metagenomics data, we have provided only 1 replicate from each sample in your instance due to limited computing and storage resources. However, if you would like to analyze complete dataset in your own computing resource, you can download all samples using the accession IDs provided in blue.

## Data location and conda environment required

Before we begin, lets remind ourselves about a few important paths and conda activation command

```bash
# Datasets

IntegratedAssignment=~/CourseData/MIC_data/IntegratedAssignment

# 16S rRNA reads
${IntegratedAssignment}/16s/16s_reads
# 16S Metadata txt
${IntegratedAssignment}/metadata.txt

# Shotgun Metagenomics
# Raw reads
${IntegratedAssignment}/ShotgunMetagenomics/metagenomic_reads

# Assembled Metagenomes
${IntegratedAssignment}/ShotgunMetagenomics/metagenomic_assemblies

# Pre-computed Metagenome Assembled Genomes (MAGs)
${IntegratedAssignment}/ShotgunMetagenomics/metagenome_assembled_genomes

# Pre-computed GTDB taxonomy assigned MAGs
${IntegratedAssignment}/ShotgunMetagenomics/GTDB_named_MAGs

# Pre-computed DRAM Annoted MAGs
${IntegratedAssignment}/ShotgunMetagenomics/DRAM_output_for_QC_MAGs/

# Activating conda environment
# Replace $software with the one you want to use
source activate $software
```

We encourage you to (not immediately) refer back to the commands provided in the presented modules, but instead to read the help functions of each function (unless you’re stuck!), you can usually bring out the help menu by using the `--help` flag.

If you really do get stuck at a step and cannot figure it out by yourself, don't worry. We left some intermediate files that you can use to "skip" a command. Simply copy it from the following locations into your working directory and continue working. Make sure you come back at a later date to troubleshoot where it went wrong!

```bash
#Intermediate files

#16s Qiime2 related files
~/CourseData/MIC_data/16s/qiime2_Intermediate_Files

#16s PiCRUSt2 related files
~/CourseData/MIC_data/16s/picrust2_Intermediate_Files
```

# PART I

## 1- 16S rRNA Gene Analyses

In this part of the integrated assignment, we will be analyzing a new 16S sequencing dataset of human oral microbiomes using what we learned in this workshop. We will be starting with FASTQ files and you will import them into Qiime2 and produce a tab-delimited file containing the ASV abundance data, as well as assign taxonomy and calculating diversity measures. The same tools will be used, but the commands to use them will not all be provided! But no worries, there will be hints throughout the assignment to point you in the right direction. One part of the exercise is to use the correct commands and parameters, and the other is to read and understand the main output files.

### 1.1 - 16S raw data preprocessing and exploration

As with all sequencing data, a sequence quality pre-check and filtering step is crucial. For your first step, let's use what we learned in the 16S module to import these FASTQ files into a Qiime2 artefact and visualize our quality scores.

After you do this, you should have a Qiime2 Artefact (*.qza) and a Qiime2 visualization file (*.qzv)

Hint commands:

```jsx
qiime tools import
qiime demux
The manifest file can be found in the 16s/16s_reads folder.
https://view.qiime2.com
```

Throughout this assignment, we will ask a few questions to keep you on track and reflect on what the outputs mean. For an added challenge, answer the questions posed in each module's practical using this oral microbiome dataset!

> **Q1. How are the sequence quality of this dataset? Are there any samples that have a low number of reads? Should the reads be truncated?**

### 1.2- ASV identification and taxonomy classification

Now that you had checked your quality, let's assign ASVs and their respective taxonomies. If you determined that you need to truncate your reads due to quality issues, remember to do that now.

At the end of this step, you should have a table of ASVs and their abundances , a table of Taxonomies assigned to each ASV, and a Qiime2 visualization artefact (*.qzv) in which you can see your Taxa bar plots.

Hint commands:

```bash
qiime dada2 denoise-paired
qiime feature-classifier
The metadata file can be found in the 16s folder.
https://view.qiime2.com
```

> **Q2. What is the difference between ASVs and OTUs? Are there major differences in taxa present at different oral sites as well as sex? Are there any site specific species?**

### 1.3- Diversity metrics & phylogenetic trees

With your ASV abundances, can you measure the diversity of our samples? Then, produce a phylogenetic tree of your microbiomes.

What you should have after this is a phylogenetic tree (*.tre) and a few Qiime2 visualization artefacts containing your diversity measurements (*.qzv)

Hint commands:

```bash
qiime diversity
qiime alignment
qiime phylogeny
https://view.qiime2.com
```

> **Q3. Using your calculated diversities, can you describe any significant between group differences?**

### 1.4- Predicting the microbiome's functional pathways using PiCRUSt2

First, in order to predict the functional pathways, we need to extract a few files from our Qiime2 workflow.

> **Q4. Do you remember what the inputs for PiCRUSt2 are?**

To generate the 3 input files from our Qiime2 work flow:

Hint commands:

```bash
#export our ASV abundances
qiime tools export --input-path unfiltered_table.qza --output-path ./
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
#removing some junk rows in the tsv file.
tail -n +2 feature-table.tsv | tail -c +2 > abundances.tsv

#export the ASV's sequences
qiime tools export --input-path representative_sequences.qza --output-path ./

```

Now that we have the ASV abundance (abundances.tsv), ASV sequences (dna-sequences.fasta) and the metadata fille (metadata.txt), let's infer some functional pathways. We are now delving into the world of PiCRUSt2, refer back to the PiCRUSt2 tutorial for a refresher if needed. Let's start by generating our predictions for ECs.

Hint commands:

```bash
place_seqs.py
hsp.py
*remember to check for NSTI values.
```

---

With the 16S and ECs predicted, lets move onto the metagenome pipeline.

Hint commands:

```bash
metagenome_pipeline.py
```

> **Q5. What do each of the outputs of this pipeline mean? Which output would you use to determine taxa shifts between functions?**

Last but not at least, pathway inference and adding descriptions.

Hint commands:

```bash
pathway_pipeline.py
add_descriptions.py
```

> **Q6. Now you have predicted the functional profile of your microbiomes, visualize it however you like and see if there are any differences in functional profiles between microbiomes of different oral sites. Can you describe the inherited limitations of using PiCRUSt2 to extrapolate microbiome functions?**

### 1.5- Conclusion

Congratulations! You have just completed the 16s portion of the integrated assignment. I hope this gave you a chance to practice and solidify your understanding of 16S marker gene metagenomics.

Inside the `~/CourseData/MIC_data/IntegratedAssignment/16s` folder, there is a file named `section1_workflow.sh` that contains all of the commands used in section 1 joined together into 1 bash script. Feel free to check your commands against it.   

Keep you findings and conclusions in mind though as we will delve into the shotgun metagenomic data of this oral microbiome in the next section. In the end, we will ask you to compare and contrast the findings of both methods.

---

---

---

# PART II

## 2- Shotgun Metegenomics

So far you have analyzed 16S part of the dataset from ***Shaiber. A et al.*** All of the samples we analyzed above were also sequenced using shotgun metagenomics sequencing as wll. The goal of this part II of assignment  is to highlight that metagenomics provide a wider lens power to view same dataset. This assignment will expand on most parts from the Module 4 and 5 tutorials.

As you have learned during the workshop, there are more than one ways to analyze metagenomics data, for this assignment we will apply both pipelines introduced in the workshop.

- **Metagenomic reads based pipeline**
- **Metagenomic assembly and binning pipeline**

## 2.1- Metagenomic reads based pipeline

### 2.1.1- Filtering out human DNA (if attached) and contamination

Since these are human oral samples in ***Shaiber. A et al.*** it's always possible that there could be some host (human) DNA that was sequenced during shotgun metagenomics. Additionally, as learned in the tutorial, sometimes PhiX sequences (a common sequencing control) aren't properly removed from output files. We'll want to remove both these types of sequences before proceeding. Process of removing the host DNA is also called as **dehosting.**

> **Q. Dehost metagenomics sample provided at** `~/CourseData/MIC_data/IntegratedAssignment/ShotgunMetagenomics/metagenomic_reads` **and prepare a summary of how many reads were dropped for each of the 6 samples?**

**Hint:** `kneaddata` from **Tutorial 4**

**Reminder:** You can use `GNU` to run samples in parallel and `--bypass-trim` to skip trimming step

### 2.1.2- Taxonomic profiling of short reads

In **Part I (1.2)** you have classified 16s rRNA ASVs using 16S based taxonomy database available. Now lets try to assign using **Whole Genome Sequencing** (WGS) based taxonomy database. Additionally, you will need to correct classification by estimating the species level abundance as discussed in Module 4 tutorial.

> **Q. Generate taxonomy profiles (classification + estimation of species level abundance) of the 6 samples that you have used in dehosting above. Do you see new taxa in metagenome based taxonomy as compared to 16S based taxonomy from Part I (1.2)? If yes, summarize if these are specie level differences/additions or wider taxa?**

**Hint:** `Kraken2` & `Bracken` from **Tutorial 4** and use the helper script provided to concatenate forward and reverse reads before running the `Kraken2`.

```bash
perl helper_scripts/concat_paired_end.pl -p 4 --no_R_match -o cat_reads kneaddata_out/*_paired_*.fastq
```

### 2.1.3- Functional profilining from short reads

Now that we have the taxonomy profiles, we want to have the potential functional profiles for these samples for better resolution.

> **Q. Map the metagenomic short reads to `UniRef90` (same database used in workshop) and summarize functional annotations across samples. Do you see unique (site-specific) functional profiles ? Are there any functional profiles that are different in both couples (B and C)?**

**Hint:** `MMseqs2` from **Tutorial 4,** you can also use `Functional_Helper_Scripts/run_TaxonomyFunctionSearchMegan.pl` to create job files for samples provided with **Tutorial 4**.

**Reminder:** As discussed in the workshop don’t forget to remove intermediate files to save significant storage space. This will help you perform following questions.

## 2.2- Metagenomic assembly and binning pipeline

This type of analyses is rather computationally intensive as compared to read based analysis. Some steps of this pipeline require significantly more CPU threads and RAM as compared to the instances available for current workshop. For this reason, we have pre-computed a few steps for you to be able to continue with next steps using the files. We have assembled metagenomes and binned Metagenome Assembled Genomes (MAGs). From these MAGs we have made 50 prominent MAGs available for further analysis. You can copy these files to your working directory and attempt to answer following questions.

### 2.2.1- Quality assesment and taxonomy classification of MAGs

Quality assesment of MAGs includes estimation of genome completeness and contamination by using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage. Assessment of genome quality can also be examined using plots depicting key genomic characteristics (e.g., GC, coding density) which highlight sequences outside the expected distributions of a typical genome.

Taxonomic classification of MAGs can be acheived by binning genomes into bins based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.

Quality assesment and taxonomic qualification are key for downstream analysis.

> **Q. Prepare a table or summary to answer following questions: Which MAGs have the highest quality in the dataset?  How many MAGs pass a >70% complete & <10% contaminated threshold?**

**Hint:** `CheckM` from **Tutorial 5**

**Reminder:** Donot forget to add `--reduced_tree` flag in CheckM, otherwise it may take significantly more time or even not possible to run `CheckM`

### 2.2.2- Prevalence, abundance, and site-specificity of MAGs

Since the MAGs were pre-computed and provided without metadata, you will need to use some information about them for further analysis. You can use your knowledge of read mapping / read recruitment from Module 5 to find out the prevalance and abundance of these MAGs in different oral sites and couples.

> **Q2. Are there any site-specific (i-e tongue-associated, plaque-associates) genomes? If so, are these genomes' species found in 16S data? Prepare a table or figure to show the prevalance and abundance of MAGs across the couples and highlight site-specific (if any).**

Hint: `Bowtie2` & `samtools` from **Tutorial 5**

### 2.2.3- Phylogenomic analysis

From first two questions in this section, we have more information about the MAGs. Now we can get into downstream analysis.

Phylogenomics (inferring genome-level evolutionary relationships) is becoming a fundamental step in many biologists’ work—such as in the characterization of newly recovered genomes, or in leveraging available reference genomes to guide evolutionary questions. You have so far built phylogenegetic tree using 16S gene in first part of this assignment, now try another one using the MAGs.

> **Q. Build a phylogenomic tree using the MAGs, can you find common taxa in phylogenetic (16S assignment Q3) and phylogenomic trees? Do you see any new taxa compared to phylogenetic tree? Are any of the MAGs potentially duplicates of each other in the tree?**

**Hint:** `GtoTree` from **Tutorial 5**

**Reminder:** For this question, we prefer for you to use MAGs that pass the QC metrics. We have already done that for you and assigned taxonomy using `GTDB`.

```bash

IntegratedAssignment=~/CourseData/MIC_data/IntegratedAssignment

# Pre-computed GTDB taxonomy assigned MAGs
${IntegratedAssignment}/ShotgunMetagenomics/GTDB_named_MAGs
```

### 2.2.5- Functional annotation of MAGs

Finally, we would like you to functionally annote QC filtered MAGs to identify different roles MAGs play within community. For functional annotation, you can use two different tool (`DRAM`, `FeGenie`). `DRAM` is computationally expensive and thus we have pre-computed the annotation process for you to visualize and summarize. Details of this are given below.

> **Q. Annotate MAGs using functional databses to identify what roles MAGs play in different oral sites? Do you see any unique functions/functional pathways in Couple B as compared to C?**

Hint: `FeGenie` & `DRAM` from **Tutorial 5**

```bash
IntegratedAssignment=~/CourseData/MIC_data/IntegratedAssignment

# Pre-computed DRAM Annoted MAGs
${IntegratedAssignment}/ShotgunMetagenomics/DRAM_output_for_QC_MAGs/
```

---

### If you are reading this part, `Congratulations` you have completed the assignment !! Lets discuss your results in our upcoming session. If you are/were stuck on any part, write us on slack `MIC` channel. If its not possible to resolve, please summarize your error for discssuion in upcoming session to discuss.
