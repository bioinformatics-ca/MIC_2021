---
layout: tutorial_page
permalink: /MIC_2021_Module6_lab
title: MIC 2021 Module 6 Lab
header1: Workshop Pages for Students
header2: Microbiome Analysis 2021
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/MIC_2021
description: MIC 2021 Module 6 Lab
---


# MetaPro Metatranscriptomics Practical Lab

**This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

This tutorial was developed by Billy Taj (billy.taj@sickkids.ca), Mobolaji Adeolu (adeolum@mcmaster.ca), John Parkinson (john.parkinson@utoronto.ca) & Xuejian Xiong (xuejian@sickkids.ca), and updated by Ana Popovic for CBW Microbiome 2021.


## Overview

This tutorial will take you through the MetaPro pipeline for processing metatranscriptomic data (https://doi.org/10.1101/2021.02.23.432558). The pipeline, developed by the Parkinson lab, consists of the following steps:

1.  Remove adapter sequences, which are added during library preparation and sequencing steps, and low-quality reads.
2.  Remove duplicate reads to reduce processing time for following steps.
3.  Remove vector contamination (reads derived from cloning vectors, spike-ins, and primers).
4.  Remove host reads (if exploring a microbiome in which the host is an issue).
5.  Remove abundant rRNA sequences which typically dominate metatranscriptomic datasets.
6.  Repopulate duplicated reads, removed in step 2.
7.  Assemble reads into contigs and predict ORFs.
8.  Annotate reads to known genes and proteins.
9.  Classify reads to known taxonomic groups and visualize the taxonomic composition of your dataset.
10.  Predict enzyme function.
11.  Generate output files, including gene annotations with normalized gene expression values, and enzyme predictions.
12.  Visualize metabolic pathway activity, as mapped onto KEGG-defined pathways, within Cytoscape.

The MetaPro metatranscriptomic pipeline includes existing bioinformatic tools and a series of Python scripts that handle the orchestrated invocation of the bioinformatics tools, read-conflict resolution, file format conversion, and output parsing. We will go through these steps to illustrate the complexity of the process and the underlying tools and scripts.  The operation of the pipeline is designed to run all steps automatically in succession.  However, for the purposes of the tutorial, a single-step tutorial-mode has been included in MetaPro's design.  One can now call the MetaPro pipeline to perform each step in isolation, for learning purposes.


New, faster, and/or more accurate tools are being developed all the time, and it is worth bearing in mind that any pipelines need to be flexible to incorporate these tools as they get adopted as standards by the community. For example, over the past two years, our lab has transitioned from Trimmomatic to AdapterRemoval, and from BLAST to DIAMOND.
Note:  This workshop was designed for use with DIAMOND v0.826.  Newer versions of DIAMOND will be incompatible with the pre-compiled database files we have made as part of this exercise.  
To illustrate the process we are going to use sequence reads generated from the contents of the colon of a mouse. These are 150 bp single-end reads. Paired-end reads can also be used, and are often preferred because they can improve annotation quality when there is enough overlap in the read pairs to improve the effective average read length. Working with paired-end data involves an additional data processing step (merging of overlapping reads) produces more files during data processing (files for merged/singleton reads, forward reads, and reverse reads), but the structure of a pipeline for paired-end data is similar to the pipeline described here and can be readily adapted.

Rather than use the entire set of 25 million read, which might take several days to process on a desktop, the tutorial will take you through processing a subset of 100,000 single-ended reads.

**Note:**
The purpose of this tutorial is to demonstrate MetaPro's various steps.  The pipeline is fully capable of running all of the steps without the intervention of the user beyond an initial call to the program.  If you wish to simply use MetaPro, do not use the --tutorial option.
This tutorial also assumes that the pipeline files are contained in a directory

## Preliminaries

### Launch the MetaPro pipeline
MetaPro operates in a containerized environment running Ubuntu Linux 18.04.  The pipeline and its dependent programs exist within a Docker container image, and may be accessed using the Docker or Singularity tools. For the purposes of this tutorial, Docker and Singlarity have both already been installed, and the MetaPro pipeline has been downloaded to your environment.  
  
Docker and Singularity maintain different access modes to use their containers.  
1) Scripted-mode: where the user calls software within the container and run. 
2) Interactive-mode: where the user can enter into the container use it like an operating systems.  

In this tutorial, we will be using Singularity in interactive mode to run MetaPro commands. The command to download and create the MetaPro image _(below)_ has already been run for you:  

&ensp;&ensp;&ensp;&ensp;singularity pull docker://parkinsonlab/metapro:develop   **[DO NOT RUN]**  


<br/><br/>Navigate to your workspace and create new tutorial folder (we will be using the absolute path):  

```
cd /media/cbwdata/workspace  
mkdir metapro_tutorial/  
cd metapro_tutorial  
```

Launch MetaPro using Singularity in interactive mode (using the command: singularity shell &lt;path to pipline metapro.sif&gt;):  

```
singularity shell /media/cbwdata/MIC_data/Module6/Tools/metapro_develop.sif
```
  
<br/><br/>To access MetaPro using Docker, we include instructions below: _(not required for the current tutorial)_  

Install Docker:  
https://www.docker.com/products/docker-desktop  

Next, pull the MetaPro docker image:  
&ensp;&ensp;&ensp;docker pull parkinsonlab/metapro:develop  

Launch MetaPro within the Docker interactive mode:  

&ensp;&ensp;&ensp;&ensp;docker run -it -v &lt;a folder in your directory>:<an equivalent folder to mount to in the container instance&gt; &lt;the docker image&gt;

&ensp;&ensp;&ensp;&ensp;An example would be:  
&ensp;&ensp;&ensp;&ensp;docker run -it -v /home/ubuntu/MetaPro_tutorial:/MetaPro_docker_tutorial parkinsonlab/metapro:develop  
  
  
  
### Download the data  

Our data set consists of 150 bp single-end Illumina reads generated from mouse colon contents. Download the data and precomputed files:  
```
wget https://github.com/ParkinsonLab/MetaPro_tutorial/releases/download/1.0/tutorial_files.tar.gz
```

Unzip the data folder, and view contents:  
```
tar -xzvf tutorial_files.tar.gz  
ls
```

The downloaded contents include:
- the input sequence file `mouse1.fastq`
- the MetaPro configuration file containing paths to tools and databases, as well as parameters `config_mouse_tutorial.ini` _(see below for description)_
- a folder containing databases required for the tutorial `databases\`
- the output directory containing some precomputed files `mouse1_run\`
- an example Cytoscape file which may be generated with MetaPro to view microbial metabolic pathway activity `Example.cys`  

MetaPro's tools may take a long time to run if the user does not have the necessary computing resources.  Therefore, we provide _pre-computed output files_ (within the `mouse1_run/` folder) so that the user is not forced to run computationally intensive steps during the tutorial. 



Change the permissions of the `mouse1_run/` folder in order to view its contents through a browser ( via your IP address, http://[ public ipv4 ] ):
```
chmod -R 777 mouse1_run
```


Access your workspace in a web browser to view the output directory: _(keep the tab open!)_    
```
http://[ insert your IPv4 ]
```



### Edit the configuration file  

MetaPro controls many of its features with a configuration file. A copy has been provided for you in the downloaded data, but it needs to be altered to include the path to the databases.

View the file:  

```
less config_mouse_tutorial.ini
```
  
Edit the configuration file to add the path to the `databases/` folder:  
- Open the file using the VI editor:  
  
```
vi config_mouse_tutorial.ini
```
- Once the file is open, navigate with arrow keys to the right of the "database_path:" field
- type `i` to activate the INSERT mode
- type the new path `/media/cbwdata/workspace/metapro_tutorial/databases`
- press ESC to exit the INSERT mode
- navigate to any characters you want to delete using arrow keys, and type `x` to remove them
- exit the file and save changes by typing a colon followed by wq `: wq` and press ENTER 
  
If you make a mistake while editing the path, you may type `u` to undo the last change, or `U` to undo all changes to the line. Alternatively, you may exit without saving any changes to the file by pressing ESC, typing `:q!` and pressing ENTER.  
  
View your modified config file:  
  
```
less config_mouse_tutorial.ini
```

```
[Databases]
database_path: /media/cbwdata/workspace/metapro_tutorial/databases
UniVec_Core: %(database_path)s/univec_core/UniVec_Core.fasta #(the UniVec core database)
Adapter: %(database_path)s/Trimmomatic_adapters/TruSeq3-PE-2.fa #(The adapters database)
Host: %(database_path)s/Mouse_cds/Mouse_cds.fasta #(The host database)
Rfam: %(database_path)s/Rfam/Rfam.cm #(The Infernal rRNA database)
DNA_DB: %(database_path)s/ChocoPhlAn/ChocoPhlAn.fasta #(The BWA database.  It assumes the index is in the same directory)
DNA_DB_Split: %(database_path)s/ChocoPhlAn/ChocoPhlAn_split/ #(The split database, for BLAT)
...
```


### Databases and licenses  

This tutorial relies on a few external databases and libraries to perform the filtering tasks associated with MetaPro.  
We have assembled the smaller databases in our precomputed files package  

[The UniVec Core database](https://ftp.ncbi.nih.gov/pub/UniVec/)  
[A mouse host sequence database](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/)  In this tutorial, we will use one from Ensembl  

There are optional databases that are mentioned in this tutorial.  Due to the size of these references, they are not required to be present for this tutorial, but if one were to use MetaPro, it is highly suggested that they are obtained:

-   [The ChocoPhlan Pangenome Database](http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/)
-   [The NCBI Non-redundant (NR) Protein Database](ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz)
-   [The GB version of the nucleotide accession2taxid table](https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)
-   [The Centrifuge NT database](https://ccb.jhu.edu/software/centrifuge/manual.shtml#nt-database)
    -   To complete the centrifuge database install, the required utilities are placed in /pipeline_tools/centrifuge
-   [The Kaiju Database](https://github.com/bioinformatics-centre/kaiju)
    -   To complete the kaiju database install, the required utilities are placed in /pipeline_tools/kaiju
    -   MetaPro relies on the full database. `(makeDB.sh -r)`
-   [The NCBI Taxdump database](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
-   [The Swiss-Prot database (fasta)](https://www.uniprot.org/downloads)
-   [The PRIAM Database](http://priam.prabi.fr/REL_JAN18/Distribution.zip)
    
These optional databases require indexing prior to use.  
-   MetaPro requires 2 versions of ChocoPhlan:
    -   One with all of the sequences in a single file, for BWA
    -   One with all of the sequences separated, for BLAT


The commands used to build the indexed databases are as follows **[DO NOT RUN]**  

    To unpack chocophlan and combine all of the sequences:
    -   tar -xvf chocophlan.tar.gz && cd chocophlan
    -   for i in $(ls | grep ".gz"); do gunzip $i; done
    -   for i in $(ls | grep ".m8"); do cat $i >> chocophlan_full.fasta
    -   mv chocophlan_full.fasta ..
     
    To prepare ChocoPhlAn for BWA:  
    -   bwa index -a bwtsw &lt;path to chocophlan_full.fasta&gt;  
    -   samtools faidx &lt;path to chocophlan_full.fasta&gt;  
                   
    To prepare NR for DIAMOND:  
    -   diamond makedb -p 8 --in &lt;path to nr&gt; -d &lt;path to nr&gt;  
                   
    To prepare Kaiju:  
    -   /pipeline_tools/kaiju/makeDB.sh -r &lt;a suitable destination for the Kaiju DB&gt;  
  
 Additionally, a license from MetaGeneMark is required to run to the contig assembly step  
 -   [MetaGeneMark](http://exon.gatech.edu/Genemark/license_download.cgi)



### Store paths to the configuration file and the output folder in variables  

```
config=/media/cbwdata/workspace/metapro_tutorial/config_mouse_tutorial.ini
output=/media/cbwdata/workspace/metapro_tutorial/mouse1_run
```

These paths will be used in running MetaPro commands which are structured as follows: `python3 /pipeline/MetaPro.py -c $config -s [sequence file] -o $output --tutorial [processing step]`  

Verify the variables:  

```
echo $config
echo $output
```

**Notes:**  

All MetaPro steps share the same file directory scheme:
- data: where the interim files are placed for each run.  This includes intermediate steps.
- final results: where the end-phase deliverables are placed, assuming the pipeline will continue running.
- All of MetaPro's commands are generated in separate shellscripts in each folder.  



### Inspect the sequences and check read quality with FastQC  

Inspect the sequences:  
```
less mouse1.fastq
```

**Notes:**
-   Type `q` to exit `less`.


Check the read quality with FastQC:  
```
/pipeline_tools/FastQC/fastqc mouse1.fastq
```

The FastQC report is generated as an HTML file `mouse1_fastqc.html`. A zip file is also generated which includes data files used to generate the report.  

Open the FastQC report in a browser!  


You can find the following information in the report:

-   Basic Statistics: Basic information of the mouse RNA-seq data, e.g. the total number of reads, read length, GC content.
-   Per base sequence quality: An overview of the range of quality values across all bases at each position.
-   Per Base Sequence Content: A plot showing nucleotide bias across sequence length.
-   Adapter Content: Provides information on the level of adapter contamination in your sequence sample.  



## Process the Reads  


### Step 1: Remove adapter sequences and trim low quality sequences.  

In the first step, MetaPro removes adaptor sequences, trims low-quality reads, and removes duplicate reads in one pass.  

The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to input sequence&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality  
  

Run the command:  

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality
```


In this Quality-filtering stage, MetaPro will perform several actions:  
- filter reads below a quality score of 75  
- filter reads below a minimum length of 30 bp  
- remove adapters  
- remove duplicate reads within the dataset  


<br/><br/>
> ***Question 1.1: How many low quality sequences have been removed?***  

<br/><br/>
Use FastQC to check the quality of the reads filtered for low quality bases and short length:  

```
cd mouse1_run/quality_filter/data/4_quality_filter/
/pipeline_tools/FastQC/fastqc singletons_hq.fastq
cd ../../../../
```

Navigate to the ~/4_quality_filter/ directory in your browser to view the HTML report.  

Compare with the previous report to see changes in the following sections:

-   Basic Statistics
-   Per base sequence quality  


<br/><br/>
> ***Question 1.2: How has the per read sequence quality curve changed in the final filtered output?***

<br/><br/>
Use FastQC to check the quality of the final filtered output:  

```
cd mouse1_run/quality_filter/final_results/
/pipeline_tools/FastQC/fastqc singletons_hq.fastq
cd ../../..
```

Compare with the previous reports to see changes in the following sections:

-   Basic Statistics
-   Per base sequence quality
-   Per sequence quality



**Notes on read quality filtering**

AdapterRemoval, which was used to remove the adapters and trim low quality bases in the reads, uses a sliding window method to remove contigous regions of low quality bases in reads. However, it is worthwhile to impose an overall read quality threshold to ensure that all reads being used in our analyses are of sufficiently error-free. For this we use the tool VSEARCH which can be found at this [website](https://github.com/torognes/vsearch) (when processing paired-end data, this step should come **after** the read merging step):

_Example command only_ **[DO NOT RUN]**. MetaPro already automatically perfoms this task.  
vsearch --fastq_filter mouse1_trim.fastq --fastq_maxee 2.0 --fastqout mouse1_qual.fastq  

The command line parameters are:
    -   `--fastq_filter ` Instructs VSEARCH to use the quality filtering algorithm to remove low quality reads
    -   `--fastq_maxee 2.0` The expected error threshold. Set at 1. Any reads with quality scores that suggest that the average expected number of errors in the read are greater than 1 will be filtered.
    -   `--fastqout` Indicates the output file contain the quality filtered reads




### Step 2. Remove duplicate reads

To significantly reduce the amount of computating time required for identification and filtering of rRNA reads, we perform a dereplication step to remove duplicated reads using the software tool CD-HIT which can be obtained from this [website](https://github.com/weizhongli/cdhit).

_Example command only_ **[DO NOT RUN]**. MetaPro already calls this command as part of the Quality-filtering step  
/usr/local/prg/cd-hit-v4.6.7-2017-0501/cd-hit-auxtools/cd-hit-dup -i mouse1_qual.fastq -o mouse1_unique.fastq


**Notes**:

The command line parameters are:
    -   `-i`: The input fasta or fastq file.
    -   `-o`: The output file containing dereplicated sequences, where a unique representative sequence is used to represent each set of sequences with multiple replicates.
A second output file `mouse1_unique.fastq.clstr` is created which shows exactly which replicated sequences are represented by each unique sequence in the dereplicated file and a third, empty, output file, `mouse1_unique.fastq2.clstr` is also created which is only used for paired-end reads.  


<br/><br/>
> ***Question 2.1: Can you find how many unique reads there are?***  

<br/><br/>Navigate to `mouse1_run/quality_filter/final_results/` to view the FastQC report, or look at the generated `singletons.fastq` file itself in the output directory.  

While the number of replicated reads in this small dataset is relatively low, with larger datasets, this step can reduce file size by as much as 50-80%  



### Step 3. Remove vector contamination

To identify and filter reads from sources of vector, adapter, linker, and primer contamination we use the Burrows Wheeler aligner (BWA) and the BLAST-like alignment tool (BLAT) to search against a database of cow sequences. As a reference database for identifying contaminating vector and adapter sequences we rely on the UniVec\_Core dataset which is a fasta file of known vectors and common sequencing adapters, linkers, and PCR Primers derived from the NCBI Univec Database. 

The format of the MetaPro command to perform vector filtering is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your quality filter final results mouse.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector  
  

Run the command: 

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/quality_filter/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector
```


MetaPro will automatically run the following:  
- Index the vector contamination database.  
&ensp;&ensp;&ensp;&ensp;bwa index -a bwtsw UniVec_Core  
&ensp;&ensp;&ensp;&ensp;samtools faidx UniVec_Core  
&ensp;&ensp;&ensp;&ensp;makeblastdb -in UniVec_Core -dbtype nucl  
- Use BWA to filter out reads aligning to the vector contamination database.  
&ensp;&ensp;&ensp;&ensp;bwa mem -t 4 UniVec_Core mouse1_unique.fastq > mouse1_univec_bwa.sam  
&ensp;&ensp;&ensp;&ensp;samtools view -bS mouse1_univec_bwa.sam > mouse1_univec_bwa.bam  
&ensp;&ensp;&ensp;&ensp;samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam  
&ensp;&ensp;&ensp;&ensp;samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam  
- Use BLAT to filter out any remaining reads that align to the vector contamination database. Sicne BLAT only accepts fasta files, VSEARCH is used to first convert the reads from fastq to fasta.  
&ensp;&ensp;&ensp;&ensp;vsearch --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta  
&ensp;&ensp;&ensp;&ensp;blat -noHead -minIdentity=90 -minScore=65  UniVec_Core mouse1_univec_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_univec.blatout  
- Lastly, a python script is used to filter the reads that BLAT does not confidently align to any vector sequences.  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/Scripts/read_BLAT_Filter_v3.py single high mouse1_univec_bwa.fastq mouse1_univec.blatout mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq


**Notes**:  

The commands do the following tasks:  
    -   `bwa index, samtools faidx, and makeblastdb`: Index the UniVec core database for BWA and BLAT  
    -   `bwa mem`: Generates alignments of reads to the vector contaminant database  
    -   `samtools view`: Converts the .sam output of bwa into .bam for the following steps  
    -   `samtools fastq`: Generates fastq outputs for all reads that mapped to the vector contaminant database (`-F 4`) and all reads that did not map to the vector contaminant database (`-f 4`)  
    
The command line parameters for BLAT are:  
    -   `-noHead`: Suppresses .psl header (so it's just a tab-separated file).  
    -   `-minIdentity`: Sets minimum sequence identity is 90%.  
    -   `-minScore`: Sets minimum score is 65. This is the matches minus the mismatches minus some sort of gap penalty.  
    -   `-fine`: For high-quality mRNAs.  
    -   `-q`: Query type is RNA sequence.  
    -   `-t`: Database type is DNA sequence.  

The argument structure for the final python script is:  
&ensp;&ensp;&ensp;&ensp;read_BLAT_Filter_v3.py &lt;operating mode: either "single" or "paired"&gt; &lt;filter stringency.  to handle paired-read conflicts.  "low" or "high"&lg; &lt;Input_Reads.fq&gt; &lt;BLAT_Output_File&lg; &lt;Unmapped_Reads_Output&gt; &lt;Mapped_Reads_Output&gt;`  

Here, BLAT does not identify any additional sequences which align to the vector contaminant database. However, we have found that BLAT is often able find alignments not identified by BWA, particularly when searching against a database consisting of whole genomes.  

In handling paired-ended data, cases will arise where one read maps to a vector, while the pair does not. The filter stringency decides how to resolve such cases:  
- Low filter stringency will only remove reads where both pairs aligned to a vector.  
- High filter stringency will remove reads where either pair aligned to a vector.  

<br/><br/>
> ***Question 3.1: Can you find how many reads BWA mapped to the vector database?***  


<br/><br/>
### Step 4. Remove host reads

To identify and filter host reads (here, reads of mouse origin) we repeat the steps above using a database of mouse DNA sequences. For our purposes we use a [mouse genome database](ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz) downloaded from Ensembl.  


The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your vector filter final results mouse.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host  
  

Run the command:  

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/vector_read_filter/final_results/singletons.fastq  
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host  
```

This call will perform the following steps:  

- Prepare the host database for alignment (BWA + BLAT)  
- perform alignment using BWA  
- convert the unaligned reads from BWA to a format for BLAT  
- perform alignment of the unaligned reads using BLAT  
- Run a script to remove the host reads from the input sample  


The following are the commands run by the script:  
&ensp;&ensp;&ensp;&ensp;bwa index -a bwtsw mouse_cds.fa  
&ensp;&ensp;&ensp;&ensp;samtools faidx mouse_cds.fa  
&ensp;&ensp;&ensp;&ensp;makeblastdb -in mouse_cds.fa -dbtype nucl  
&ensp;&ensp;&ensp;&ensp;bwa mem -t 4 mouse_cds.fa mouse1_univec_blat.fastq > mouse1_mouse_bwa.sam  
&ensp;&ensp;&ensp;&ensp;samtools view -bS mouse1_mouse_bwa.sam > mouse1_mouse_bwa.bam  
&ensp;&ensp;&ensp;&ensp;samtools fastq -n -F 4 -0 mouse1_mouse_bwa_contaminats.fastq mouse1_mouse_bwa.bam  
&ensp;&ensp;&ensp;&ensp;samtools fastq -n -f 4 -0 mouse1_mouse_bwa.fastq mouse1_mouse_bwa.bam  
&ensp;&ensp;&ensp;&ensp;vsearch --fastq_filter mouse1_mouse_bwa.fastq --fastaout mouse1_mouse_bwa.fasta  
&ensp;&ensp;&ensp;&ensp;blat -noHead -minIdentity=90 -minScore=65  mouse_cds.fa mouse1_mouse_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_mouse.blatout  
&ensp;&ensp;&ensp;&ensp;./1_BLAT_Filter.py mouse1_mouse_bwa.fastq mouse1_mouse.blatout mouse1_mouse_blat.fastq mouse1_mouse_blat_contaminats.fastq  


<br/><br/>
> ***Question 4.1: How many reads did BWA and BLAT align to the mouse host sequence database?***  



***Optional:*** In your own future analyses you can choose to complete steps 3 and 4 simultaneously by combining the vector contamination database and the host sequence database using `cat UniVec_Core mouse_cds.fa > contaminants.fa`. However, doing these steps together makes it difficult to tell how much of your reads came specifically from your host organism.  


### Step 5. Remove abundant rRNA sequences *** **[DO NOT RUN]**  

rRNA genes tend to be highly expressed in all samples and must therefore be screened out to avoid lengthy downstream processing times for the assembly and annotation steps. MetaPro uses [Barrnap](https://github.com/tseemann/barrnap) and [Infernal](http://infernal.janelia.org/).
You could use sequence similarity tools such as BWA or BLAST for this step, but we find Infernal, albeit slower, is more sensitive as it relies on a database of covariance models (CMs) describing rRNA sequence profiles based on the Rfam database. Due to the reliance on CMs, Infernal, can take as much as 4 hours for ~100,000 reads on a single core.  In an effort to shrink the computing time, we leverage a computing cluster's multiple cores.  

Here, MetaPro demonstrates the case for automation. MetaPro subdivides the input data, coordinates the concurrent processes, and collects the results into one single file after all of the scanning has been complete.  


MetaPro will perform the following:  

- subdivide the input data into user-defined chunk sizes (e.g. 1000 reads).  
- Each chunk is then run independently:  
  - run each chunk through Barrnap  
  - using the results of Barrnap, filter the data chunk into mRNA, and leftover data for further scanning.  
  - run each leftover chunk through Infernal.  
  - filter the Barrnap leftover chunk using the Infernal results, to get mRNA, and "other"  
- collect all of the data pieces (Barrnap mRNA, Infernal mRNA) into mRNA, and "other"  


By running things this way, the rRNA step takes 4 minutes (as recorded with a 40-core computing node with 200 GB RAM, and an rRNA chunksize of 1000 reads), but it requires significant computing power, memory, and storage space, not available on a typical desktop PC.  

If you were to run this on your own, you will need the RFam database.  

The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your host filter final results mouse.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA  


The command would look as follows: **[DO NOT RUN]**  

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/host_read_filter/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA
```

We have provided you with the pre-computed results:  

``` 
ls mouse1_run/rRNA_filter/final_results
```

Here, we only remove a few thousand reads that map to rRNA, but in some datasets rRNA may represent up to 80% of the sequenced reads.  

<br/><br/>
> ***Question 5.1: How many rRNA sequences were identified? How many reads are now remaining?***  

<br/><br/>
### Step 6. Rereplication / duplicate repopulation

After removing contaminants, host sequences, and rRNA, we need to replace the previously removed replicate reads back in our data set.  

The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your rRNA filter final results mouse.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop  


Run the command:  

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/rRNA_filter/final_results/mRNA/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop
```

Now that we have filtered vectors, adapters, linkers, primers, host sequences, and rRNA, check read quality of putative mRNA reads with FastQC:  

```
/pipeline_tools/FastQC/fastqc mouse1_run/duplicate_repopulation/final_results/singletons.fastq
```  

<br/><br/>
> ***Question 6.1: How many total contaminant, host, and rRNA reads were filtered out?***  

<br/><br/>
### Step 7. Contig assembly   *** **[DO NOT RUN]**

We have now identified the putative mRNA reads, and here we assemble the reads into contigs. Previous studies have shown that assembling reads into larger contigs significantly increases the ability to annotate them to known genes through sequence similarity searches. Here we will apply the SPAdes transcript assembly algorithm to our set of putative mRNA reads.  

The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your rereplicated mRNA final results mouse.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs  
  

The command would look as follows: **[DO NOT RUN]**  
  
```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/duplicate_repopulation/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs
```

The pre-computed output is located in the following directory:  
```
ls mouse1_run/assemble_contigs/final_results
```


**Notes**:
In this step, MetaPro does the following:  

-   Assemble the reads into contigs using SPAdes
-   Use MetaGeneMark to predict the genes in these contigs
-   Use BWA to align the mRNA reads against these split-contigs to find out which reads were consumed by the process.
-   Produce a relational map of split-contig and their constituent reads.

[SPAdes](https://cab.spbu.ru/software/spades/) assembles long contigs, but MetaPro requires that each contig only represent 1 gene.  Thus the need to disassemble them, using [MetaGeneMark](http://exon.gatech.edu/Genemark/meta_gmhmmp.cgi)  MetaGeneMark requires the user to register and obtain a free license.  Thus, we have provided the results in the precomputed files package.  


<br/><br/>
> ***Question 7.1: How many contigs did SPAdes produce?  
Hint: try using the command `tail mouse1_contigs.fasta`***  
  
<br/><br/>
> ***Question 7.2: How many reads were not used in contig assembly? How many reads were used in contig assembly?***  
  

<br/><br/>
### Step 8. Annotate reads to known genes/proteins *** **[DO NOT RUN]**  


Here we will attempt to infer the specific genes our putative mRNA reads originated from. In our pipeline we rely on a tiered set of sequence similarity searches of decreasing accuracy - BWA, BLAT, and DIAMOND. While BWA provides high stringency, sequence diversity that occurs at the nucleotide level results in few matches observed for these processes. Nonetheless it is quick. To avoid the problems of diversity that occur at the level of nucleotide, particularly in the absence of reference microbial genomes, we use a cascaded method involving two other tools: BLAT, and DIAMOND. BLAT provides a more sensitive alignment, along with quality scores to rank the matches.  DIAMOND is used to provide more sensitive peptide-based searches, which are less prone to sequence changes between strains.

Since BWA and BLAT utilize nucleotide searches, we rely on the [ChocoPhlan pangenome database](http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/) that we obtained from The Huttenhower lab, which contains over 10000 organisms in separate .ffn files.  We create a merged copy of these sequences, and index it for BWA to use.  We leave it in its separated state for BLAT to use.  

For DIAMOND searches we use the [Non-Redundant (NR) protein database](ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz) from the NCBI.

This is a computationally intensive step.  We employ our subdivision strategy here, similar to our design for rRNA removal, as seen below  

-   The mRNA read data (contigs, remaining singletons, and remaining paired reads if applicable) are split into chunksizes (GA_chunksize in the configuration)
-   Each chunk is sent through BWA to be aligned against the ChocoPhlan database  
-   A map of genes to constituent reads is formed for each chunk.
-   Reads not annotated by BWA are isolated, and are sent through to BLAT.
-   Another map of genes to constituent BLAT-annotated reads are formed for each chunk.
-   Reads not annotated by BLAT are isolated, and are sent through to DIAMOND.
-   A 3rd batch of maps of annotations to constituent reads are formed for the DIAMOND reads for each chunk.

In all, this leaves us with (split into chunks):  

-   a batch of BWA-annotated gene-to-read maps
-   a batch of BLAT-annotated gene-to-read maps
-   a batch of DIAMOND-annotated protein-to-read maps
-   a batch of reads unannotated by BWA, BLAT, and DIAMOND.

These batches of files are then sent to a custom script that will perform all of the final merging:  

-   collect and merge every map to form one map.  
-   collect and merge the unannotated DIAMOND reads
-   collect all of the genes that were found in the reads, and convert them into proteins.  Then merge them with the proteins found in DIAMOND.  This step is for downstream analysis.  


The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your unassembled singletons.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;contig='&lt;path to your contigs.fasta&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA  
  
  
The command would look as follows: **[DO NOT RUN]**  _This step is heavily computationally intensive._  
  
```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA
```

We have provided pre-computed outputs in the following directory:  
  
```
ls mouse1_run/GA_FINAL_MERGE/final_results
```


**Notes**:
-   The gene annotation step will create 4 different subdirectories: GA_BWA, GA_BLAT, GA_DIAMOND, and GA_FINAL_MERGE
-   MetaPro will take the top-quality hits of each tool to count towards annotation:
    -   In BWA: the CIGAR string is decoded.  Any hit with less than 90% match is rejected 
    -   In BLAT and DIAMOND: the only reads that pass annotation are reads with all 3 conditions satisfied:
        -   A sequence identity score of greater than 85
        -   An alignment length of greater than 65%
        -   A bitscore of over 60
-   MetaPro makes a number of extra considerations to account for multi-mapped reads:
    -   Every read scanned by BWA, BLAT, and DIAMOND is affixed with a quality score suffix taken from the aligner's report
    -   In BWA: the alignment score is used.
    -   In BLAT and DIAMOND: the bitscore is used.
    -   Should there be a case where a read is annotated to multiple genes or proteins, each tool will use their quality scores to rank the hit
        -   Alignment information is iterated through, from the top of the file down to the bottom.
        -   The hit with the higher score is used for all cases.
        -   If the scores are tied, then the incumbent hit is used.
    -   This read-ranking is done on each chunk when the gene-to-read maps are formed.
    -   There are further considerations made for paired-end annotation:
        -   MetaPro views paired-end data as 2 copies of the same read
        -   If there are disagreements between the forward and reverse read's annotation, the quality scores are used to resolve the conflict.
        -   If both reads agree on the same gene or protein, the read is counted only once. 
        -   This paired-read conflict resolution is performed in the GA_FINAL_MERGE step.
    
-   Unless you are running this tutorial on a computing cluster, most systems do not have enough memory to handle indexing or searching large databases like `ChocoPhlan` (19GB) and `nr` (>60GB). The descriptions in this section are purely for your information. Please use our precomputed gene, protein, and read mapping files from the tar file `tar -xzf tutorial_files.tar.gz`  
-   

### Step 9. Taxonomic Classification *** **[DO NOT RUN]**  

Now that we have putative mRNA transcripts, we can begin to infer the origins of our mRNA reads. Firstly, we will attempt to use a reference based short read classifier to infer the taxonomic orgin of our reads. Here we will use [Kaiju](https://github.com/bioinformatics-centre/kaiju), [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml), and our Gene Annotation results to generate taxonomic classifications for our reads based on a reference database. 
Kaiju can classify prokaryotic reads at speeds of millions of reads per minute using the proGenomes database on a system with less than 16GB of RAM (~13GB). Using the entire NCBI nr database as a reference takes ~43GB. Similarly fast classification tools require >100GB of RAM to classify reads against large databases. 

However, Kaiju still takes too much memory for the systems in the workshop so we have precompiled the classifications, `mouse1_classification.tsv`, in the tar file `tutorial_files.tar.gz`.
Centrifuge is a lightweight rapid microbial classification engine.  It uses methods similar to BWA and the Ferrgina-Manzini (FM) index  to make quick work of assigning taxomony.

The ChocoPhlan Pangenome Database contains taxonomic information that MetaPro extracts.  Kaiju, Centrifuge, and the extracted taxa are combined using [WEVOTE](https://github.com/aametwally/WEVOTE).  WEVOTE is the Weighted Voting Taxonomic Identification system.  It performs consensus merging of various taxa results and reconciles the taxa identification from various sources.  
MetaPro uses this to settle on one confident taxon amongst Kaiju, Centrifuge, and the ChocoPhlan database choices.  


The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your unassembled singletons.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;contig='&lt;path to your contigs.fasta&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA  
  

The command would look as follows: **[DO NOT RUN]**  _MetaPro assumes the Gene Annotation step has completed._  
  
```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA
```


We have provided pre-computed results here:
```
mouse1_run/taxonomic_annotation/final_results
```

We can use [Krona](https://github.com/marbl/Krona/wiki) to generate a hierarchical multi-layered pie chart summary of the taxonomic composition of our dataset.  First, the export of MetaPro's taxonomic annotations needs to be slightly modified. Run the following commands to generate a Krona output:  

```
python3 /pipeline/Scripts/alter_taxa_for_krona.py mouse1_run/taxonomic_annotation/final_results/taxonomic_classifications.tsv mouse1_classification.tsv
/pipeline_tools/kaiju/kaiju2krona -t databases/nodes.dmp -n databases/names.dmp -i mouse1_classification.tsv -o mouse1_classification_Krona.txt
/pipeline_tools/KronaTools/scripts/ImportText.pl -o mouse1_classification.html mouse1_classification_Krona.txt
```

View the pie chart representation of the taxonomies detected through a web browser.  

<br/><br/>
> ***Question 9.1: What is the most abundant family in our dataset? What is the most abundant phylum?  
Hint: Try decreasing the `Max depth` value on the top left of the screen and/or double clicking on spcific taxa.***


<br/><br/>
### Step 10. Enzyme Function Annotation *** **[DO NOT RUN]**

To help interpret our metatranscriptomic datasets from a functional perspective, we rely on mapping our data to functional networks such as metabolic pathways and maps of protein complexes. Here we predict enzyme functions (EC numbers) in our dataset, such that we may later map them to KEGG metabolic pathways.

MetaPro uses three tools to produce its enzyme annotations: [DETECT](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq266), [PRIAM](http://priam.prabi.fr/), and DIAMOND.


We use DIAMOND to identify homologs of our genes/proteins in the SWISS-PROT database that have assigned enzyme functions. This is a relatively coarse and straightforward way to annotate enzyme function by homology. We also use more robust methods for enzymatic function annotation, such as our own probability density-based enzyme function annotation tool, DETECT, and the tool PRIAM. MetaPro combines the predictions of all three tools to give two answers, a lower-confidence set of enzyme predictions, and a higher-confidence set of predictions.  

PRIAM is incredibly resource-intensive and slow to run.  For the sake of time, pre-computed results have been provided.


The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your unassembled singletons.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;contig='&lt;path to your contigs.fasta&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC  


The command would look as follows: **[DO NOT RUN]** _This command assumes that MetaPro has performed the Gene annotation steps._  

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC
```

The pre-computed results are provided in the following directory:  
```
ls mouse1_run/enzyme_annotation/final_results
```



**Notes:**  
MetaPro's high-confidence and low-confidence are determined by the following:  
    -   DIAMOND: Low-confidence hits are ones with an e-value of 1e-5 or smaller.  High-confidence hits are ones with an e-value of 1e-10 or smaller  
    -   PRIAM: Low-confidence hits are ones with e-values lower than 1e-5.  High-confidence hits are ones where their probability value is 0.5 or higher  
    -   DETECT: There is no separation.  
    
MetaPro reconciles the 3 annotations in the following manner:  
    -   The Enzymes predicted by DETECT are taken, followed by the annotations that agree between PRIAM and DIAMOND
    -   In the event of multiple enzymes being annotated to the same protein:
        -   Every enzyme annotation comes with a probability score.
        -   MetaPro includes a enzyme co-occurence database (compiled from the ENZYME database) that contains pairs of enzymes known to exist together
        -   Using this co-occurence database, MetaPro filters invalid predictions.
            -   In cases where more-than-2 enzymes are annotated to a protein, the top 2 enzymes are taken, based on the probability score.
            -   If a pair of enzymes do not exist in the database, the enzyme with the higher probability score is declared the proper annotation.
            -   Otherwise, the annotation is declared as a pair of enzymes.  
  

<br/><br/>
> ***Question 10.1: How many high-confidence unique enzyme functions were identified in our dataset?***  


<br/><br/>
### Step 11. Generate output files *** **[DO NOT RUN]**  

We have removed low quality bases/reads, vectors, adapters, linkers, primers, host sequences, and rRNA sequences and annotated reads to the best of our ability. We will now generate summaries of the gene counts, predicted functions and taxonomies in our microbiome.  

MetaPro generates many output files:  

-   An account of read numbers during various filtering and annotation steps
-   A gene expression table of counts and RPKM values 
-   RPKM values of genes in the 20-most prevalent taxa in the sample  
-   A Cytoscape-compatible network file  
-   An enzyme superpathway heatmap to visualize the distribution of enzymes found  
-   A gene/protein-to-read map of all genes and proteins identified by MetaPro, followed by its constituent reads  
-   A histogram of read quality  
-   A summary of all taxa identified, followed by the number of reads associated with those taxa  


The format of the MetaPro command is:  

&ensp;&ensp;&ensp;&ensp;read1='&lt;path to your unassembled singletons.fastq&gt;'  
&ensp;&ensp;&ensp;&ensp;contig='&lt;path to your contigs.fasta&gt;'  
&ensp;&ensp;&ensp;&ensp;python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial output  


The command would appear as follows:  **[DO NOT RUN]**  _The command assumes that MetaPro has performed the gene, taxa, and enzyme annotations._  

```
read1=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/cbwdata/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial output  
```

We have provided the pre-computed output files generated by MetaPro. View the files:  

```
ls mouse1_run/outputs/final_results/
```

<br/><br/>
> ***Question 11.1: How many unique enzyme activities were predicted, at low and high confidence? View the `read_count.txt` file.***  

<br/><br/>
> ***Question 11.2: Have a look at the `RPKM_table.tsv` file. What are the most highly expressed genes?***  

<br/><br/>
> ***Question 11.3: Which taxa appear to contribute the most metabolic activity? View the `enzyme_superpathway_heatmap.jpg`.***  


<br/><br/>
### Step 12. Visualize the results using a KEGG Pathway as a scaffold in Cytoscape.  

This step will be performed on **your workstation**, using a Cytoscape file output from MetaPro.  

To visualize our processed microbiome dataset in the context of the carbohydrate metabolism pathways, we use the network visualization tool **Cytoscape** together with the `enhancedGraphics` and `KEGGscape` plugins. Some useful commands for loading in networks, node attributes and changing visual properties are provided below (there are many Cytoscape tutorials available online).  Use [Cytoscape 3.7.2](https://github.com/cytoscape/cytoscape/releases/3.7.2/) instead of the latest version.  


**Download the metabolic pathway**

First, download the carbohydrate metabolism pathways from KEGG onto your workstation by pasting the following addresses into a browser:  
```
https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml
```


Alternatively, if you are using Linux, you may use the `wget` command as follows:
```
wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml
```

You can find other [pathways on KEGG](http://www.genome.jp/kegg-bin/get_htext?htext=br08901.keg) which can also be imported into Cytoscape by selecting the `Download KGML` option on the top of the page for each pathway.


Next, download the `Cytoscape_network.tsv` file through the browser. It is located within the `mouse1_run/outputs/final_results` folder.  This file contains expression values of enzymes (as RPKMs) of all taxa detected at >1% abundance within the dataset.


**Install the Cytoscape plugins**

-   Select `Apps` -> `App Manager`
-   Search for `enhancedGraphics`
-   Select `enhancedGraphics` in the middle column then click `Install` in the bottom right
-   Search for `KEGGScape`
-   Select `KEGGScape` in the middle column then click `Install` in the bottom right

**Import an XML from KEGG into Cytoscape**

-   Select `File` -> `Import` -> `Network` -> `File...`
-   Select the XML file, `ec00010.xml` or `ec00500.xml` and click `Open`
-   Check `Import pathway details from KEGG Database` box then select `OK`

**Loading a node attribute text file (.txt) - this will map attributes to nodes in your network which you can subsequently visualize**

-   Select `File` -> `Import` -> `Table` -> `File...`
-   Select the `Cytoscape_network.tsv` file and click `Open`
-   Change the `Key Column for network` from `shared name` to `KEGG_NODE_LABEL`
-   Click OK

**Visualizing your node attributes**

-   In the left `Control Panel` select the `Style` tab
-   Check the `Lock node width and height` box
-   Click the left-most box by the `Size` panel and change the default node size to 20.0
-   Click the blank box immediately to the right of the box you clicked to change the default size, change the `Column` field to `RPKM` and the `Mapping Type` field to `Continuous Mapping`
-   Click the left-most box by the `Image/Chart 1` panel, switch to the `Charts` tab, Click the doughnut ring icon, and press the `>>` "add all" button between the two column fields before clicking apply (make sure to remove overall RPKM from the fields that are added to the doughnut ring)
-   If you do not see the `Image/Chart 1` panel, select `Properties` -> `Paint` -> `Custom Paint 1` -> `Image/Chart 1` from the to left corner of the control panel
-   To improve the visualization you can modify colour properties under `Image/Chart 1` -> `Charts` -> `Options`, or modify other properties such as Label Font Size, Label Position, Fill Color, Node location, and edge properties

**Notes:**

-   A cytoscape file with node attributes precalculated is provided for your convenience, `tar -xzf tutorial_files.tar.gz Example.cys`, feel free to open it and play with different visualizations and different layouts - compare the circular layouts with the spring embedded layouts for example. If you want to go back to the original layout then you will have to reload the file.
-   Cytoscape can be temperamental. If you don't see pie charts for the nodes, they appear as blank circles, you can show these manually. Under the 'properties' panel on the left, there is an entry labeled 'Custom Graphics 1'. Double click the empty box on the left (this is for default behavior) - this will pop up a new window with a choice of 'Images' 'Charts' and 'Gradients' - select 'Charts', choose the chart type you want (pie chart or donut for example) and select the different bacterial taxa by moving them from "Available Columns" to "Selected Columns". Finally click on 'Apply' in bottom right of window (may not be visible until you move the window).

**Visualization Questions:**

- Which genes are most highly expressed in these two systems?
- Which taxa are responsible for most gene expression?
- Can you identify sub-systems (groups of interacting genes) that display anomalous taxonomic profiles?
- Think about how you might interpret these findings; for example are certain taxa responsible for a specific set of genes that operate together to fulfill a key function?
- Can you use the gene annotations to identify the functions of these genes through online searches?
- Think about the implications of sequence homology searches, what may be some caveats associated with interpreting these datasets?
