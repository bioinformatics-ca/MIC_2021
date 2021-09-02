---
layout: tutorial_page
permalink: /MIC_2021_Module4_lab
title: MIC 2021 Module 4 Lab
header1: Workshop Pages for Students
header2: Microbiome Analysis 2021
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/MIC_2021
description: MIC 2021 Module 4 Lab
---
This tutorial is for the 2021 Canadian Bioinformatics Workshop running from Sept. 1st -Sept. 3rd.

**Authors**: Jacob Nearing and Morgan Langille

This tutorial has been adapted from the previous 2018 version written by Morgan Langille and Gavin Douglas.

### Table of Contents

## Introduction

The main goal of this tutorial is to introduce students to different approaches for the taxonomic profiling and functional profiling or metagenomic data from microbiome samples. Throughout this tutorial we will emphasize that there is not a one-size-fits-all pipeline for analyzing MGS data. For example some individuals may opt to examine their data using a marker-based approach such as MetaPhlAn3 while others may opt to use a kmer based strategy such as Kraken2+Braken. Furthermore, in a subsequent module in this workshop students will also learn about another approach for examine microbiome data using metagenomic assembly and binning.

Throughout this tutorial there will be questions are aimed at helping students understand the various steps in reference based metagenomic profiling. [The answers are found on this page]

### Teaching Objectives
1. Become familiar with raw shotgun metagenomics data
2. Learn how to use [GNU Parallel](https://www.gnu.org/software/parallel/) to easily run commands on multiple files
3. Learn common approaches to pre-processing MGS data
4. Learn how to generate taxonomic profiles from your MGS data using both kmer based strategies [Kraken2]() + [Braken]().
5. Learn how to generate functional profiles using MMseqs2 in combination with Kraken2.

### Bioinformatic tool citations
Properly citing bioinformatic tools is important to ensure that developers are getting the proper recognition. If you use any of the commands in this tutorial for your own work be sure to cite the relevant tools listed below!

* **Bowtie2** ([website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/))
* **FASTQC** ([website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
* **GNU Parallel** ([website](https://www.gnu.org/software/parallel/), [paper](https://www.usenix.org/system/files/login/articles/105438-Tange.pdf))
* **kneaddata** ([website](https://bitbucket.org/biobakery/kneaddata/wiki/Home))
* **Microbiome Helper** ([website](https://github.com/LangilleLab/microbiome_helper/wiki), [paper](http://msystems.asm.org/content/2/1/e00127-16))
* **MinPath** ([website](http://omics.informatics.indiana.edu/MinPath/), [paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000465))
* **Kraken2** ([website](https://ccb.jhu.edu/software/kraken2/), [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)
* **Bracken** ([website](https://ccb.jhu.edu/software/bracken/), [paper](https://peerj.com/articles/cs-104/)
* **MMseqs2** ([website](https://github.com/soedinglab/mmseqs2), [paper](https://www.nature.com/articles/nbt.3988)

### Activate conda environment
[Anaconda] or for short conda is a power tool used to managed different environments within your computer. Many bioinformatic tools can be installed using conda as it can significantly reduce the amount of leg work required when setting up the proper dependencies and versioning of tools. As such in this tutorial we have included a conda environment that contains all of the required bioinformatic tools for this tutorial!

You can activate this environment using the following command:

```
source activate CBW_2021_KMER
```
When you are finished this tutorial you can deactivate the conda environment using:

```
conda deactivate
```

### Introducing GNU Parallel
Often in bioinformatics research we need to run files corresponding to different samples through the same pipeline. You might be thinking that you can easily copy and paste commands and change the filenames by hand, but this is problematic for two reasons. Firstly, if you have a large sample size this approach will simply be unfeasible. Secondly, sooner or later you will make a typo! Fortunately, [GNU Parallel](https://www.gnu.org/software/parallel/) provides a straight-forward way to loop over multiple files. This tool provides an enormous number of options: some of the key features will be highlighted below.

The first step will be to create a directory that we will be working in.

```
cd workspace
mkdir Module4
cd Module4
```

Next we will need to download the zip folder containing our test data using the command ```wget```:

```
wget https://www.dropbox.com/s/58grzx1ir7o8d3k/gnu_parallel_examples.zip?dl=1 -O gnu_parallel_examples.zip
unzip gnu_parallel_examples.zip
``` 

We can now safely remove the zip file and move into the gnu_parallel_examples folder.
```
rm gnu_parallel_examples.zip
cd gnu_parallel_examples
```

Firstly, whenever you're using ```parallel``` it's a good idea to use the ```--dry-run``` option, which will print the commands to screen, **but will not run them**.

This command is a basic example of how to use ```parallel```:

```
parallel --dry-run -j 2 'gzip {}' ::: example_files/fastqs/*fastq
```

This command will run ```gzip``` (compression) on the test files that are in ```example_files/fastqs/```. There are 3 parts to this command:

1. The options passed to GNU Parallel: ```--dry-run -j 2```.
    The ```-j``` option specifies how many jobs should be run at a time, which is 2 in this case. This means that both of the commands printed out will be run at the same time. This is what makes parallel so powerful as it allows you to run multiple commands on different files at the same time. 

**NOTE you should avoid trying to run more jobs at once than there are available processors on your workstation**
The workstations that we are working on have a total of 4 CPUs! As such we should not run parallel commands with more than 4 jobs at a time. This will vary depending on the hardware you are working on. 

2. The commands to be run: ```'gzip {}'```. The syntax ```{}``` stands for the input filename that is the last part of the command. In parallel the command you want to run should always be surrounded by the character ```'```

3. The files to loop over: ```::: example_files/fastqs/*fastq```. Everything after ```:::``` is interpreted as the files to read in!

Remove ```--dry-run``` from the above command and try running it again. This should result in both fastq files located in ```example_files/fastqs/``` to be compressed in gzip format. You can recognise this compression type as the files will always end in ```.gz```.

Sometimes multiple files need to be in given to the same command you're trying to run in parallel. As an example of how to do this, take a look at the paired-end FASTQ files in this example folder:

```
ls example_files/paired_fastqs
```

The "R1" and "R2" in the filenames specify forward and reverse reads respectively. If we wanted to concatenate the forward and reverse files together (paste them one after another into a single file for each pair) we could use the below command:
```
parallel --dry-run -j 2 --link 'cat {1} {2} > {1/.}_R2_cat.fastq' ::: example_files/paired_fastqs/*R1.fastq ::: example_files/paired_fastqs/*R2.fastq
```

There are a few additions to this command:
1. Two different sets of input files are given, which are separated by ```:::``` (the first set is files matching *R1.fastq and the second set is files matching *R2.fastq).
2. The ```--link``` option will allow files from multiple input sets to be used in the same command in **order**. ```{1}``` and ```{2}``` stand for a input file from the first and second set respectively.
3. Finally, the output file was specified as ```{1/.}_R2_cat.fastq```. The ```/``` will remove the full path to the file and ```.``` removes the extension.

It's important to check that each samples' files are being linked correctly when the commands are printed to screen with the ```--dry-run``` option. If the commands look ok then try running them without ```--dry-run``` now!

You should now have two files in your directory named ```test[1-2]_R1_R2_cat.fastq```. Lets see how many lines these files contain with the following command:

```
wc -l test*
```

You can return to the higher directory with this command:

```
cd ..
```

Finally one thing to keep in mind when running slow command with ```parallel``` is that you can keep track of their progress (e.g. how many are left to run and the average runtime) with the ```--eta``` option.

### Exploring the raw data

Alright now let's take a look at the data we'll be using for the main tutorial! The input files are FASTQ files, which contain information on the base-pairs in sequenced reads and the quality scores of each base-pair.

One thing we should be aware of when dealing with metagenomic data is that many tools will require a large amount of computational resources in order to run. As such they can either take awhile to run or in some cases fail to run at all due to issues such as insufficient memory. In many cases it is unlikely that metagenomic analysis of microbiome samples can be completed on a standard laptop in a timely manner.  

For now we'll work with two example sub-sampled FASTQs, that were downloaded from the [HMP2 database that examined individuals with and without inflammatory bowel disease](https://www.hmpdacc.org/ihmp/). These files are much smaller than normal metagenomic sequencing files but will speed up our data processing during this tutorial!

First we will need to copy the files from the storage section of our machine into our workspace.

```
cp -r ~/CourseData/MIC_data/Module4/Data/raw_data/ .
```
Note the -r flag on the copy command indicates that we want to copy a directory and not an individual file.



Take a look with ```less``` at the forward ("R1") and reverse ("R2") reads. 

```
less raw_data/CSM79HR8_R1_subsampled.fastq.gz
less raw_data/CSM79HR8_R2_subsampled.fastq.gz
```

You'll see that there are 4 lines for each individual read. Also, note that in the forward and reverse FASTQs for each sample that the read ids (the 1st line for each read) are the same and in the same order (except for ```/1``` or ```/2``` at the end).

To find out how many reads are in each forward FASTQ we can simply count the number of lines in a file and divide it by 4. Since the files are gzipped we'll need to use ```zcat``` to uncompress them on the fly and pipe (```|```) this as input to the ```wc -l``` command, which will return the number of lines in plaintext files.

```
zcat raw_data/CSM79HR8_R1_subsampled.fastq.gz | wc -l
zcat raw_data/HSM7J4QT_R1_subsampled.fastq.gz | wc -l
```
We can see that the files have a similar number of lines and therefore have very similar read depths. However, this is not always the case with MGS data!
 
**Question 1: Based on the above line counts can you tell how many reads should be in the reverse FASTQs? Remeber that a standard convention is to mark the reverse FASTQs with the characters ```R2```**

## Pre-processing

Pre-processing is an important step when analyzing sequencing data. When not performed you could be allowing contaminants (such as DNA from the host) and/or low-quality reads to introduce noise to downstream analyses.

### Visualizing read quality

A standard first step with pre-processing sequence data is to visualize the base-qualities over the reads along with other helpful statistics. FASTQC is a simple program that will allow you to explore this information!

You can run FASTQC on the two test samples with this command (after making the output directory):

```
mkdir fastqc_out
fastqc -t 4 raw_data/*fastq.gz -o fastqc_out
```
This is a fairly simple tool and only has a few arguements.
* ```-t``` - The number of threads we want to use to process each file
* ```raw_data/*fastq.gz``` - This indicates we want to run this command on all files that end with ```fastq.gz``` in the folder raw_data
* ```-o``` - Indicates the name of the folder we want to output the results to. 

You can open one of the HTML files in your browser by entering the same address you used to connect to your AWS instance into your favourite web browsers. You should see something like this under the per-base quality section:

[[images/CBW_2018/CSM79HR8_R1_fastqc_qual.png]]

These base qualities are **not** typical! Usually you would see the quality scores drop off near the end of the reads (and this tends to be even worse for reverse reads). However, it turns out that these reads have already been filtered! Importantly, this wouldn't have been obvious without visualising the quality scores.

### Filtering out low-quality reads

Even though these reads have already been through one round of quality filtering, we'll check them for contaminant sequences. Since these are human stool samples it's possible that there could be some human DNA that was sequenced (and could be high quality!). Also, sometimes PhiX sequences (a common sequencing control) aren't properly removed from output files. We'll want to remove both these types of sequences before proceeding.

We'll use ```kneaddata```, which is a wrapper script for multiple pre-processing tools, to check for contaminant sequences. Contaminant sequences are identified by this tool by mapping reads to the human and PhiX reference genomes with Bowtie2. You would also typically trim trailing low-quality bases from the reads at this step with Trimmomatic, which kneaddata will do automatically if you didn't use the ```skip-trim``` option below.

The below command will pass the read pairs to the same command, similar to how we ran them above. If you want to run quality processing as well checkout the command in our [metagenomics SOP](https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Standard-Operating-Procedure-v3). Note that in this case we're running 1 job at a time (```-j 1```) and each job is being run over 4 threads (```-t 4```). This means we are using 4 processors (CPUs) at a time while running this program (the maximum number on your workstation!) The Bowtie2 index files have already been made and are passed to the tool with this option: ```-db ~/CourseData/MIC_data/Module4/Data/bowtie2_db/GRCh38_PhiX```.


```
parallel --eta -j 1 --link 'kneaddata -i {1} -i {2} -o kneaddata_out/ \
-db  /home/ubuntu/CourseData/MIC_data/Module4/Data/bowtie2db/GRCh38_PhiX \
-t 4 --bypass-trim ' \
 ::: raw_data/*_R1_subsampled.fastq.gz ::: raw_data/*_R2_subsampled.fastq.gz
```

If you take a peek in ```kneaddata_out``` you'll see that a number of output files have been created! The key information we want to know is how many reads were removed. This command will generate a simple count table:

```
kneaddata_read_count_table --input kneaddata_out --output kneaddata_read_counts.txt
```

**Question 2: Take a look at ```kneaddata_read_counts.txt```: how many reads in sample CSM79HR8 were dropped?**

You should have found that only a small number of reads were discarded at this step, which is **not** typical. It happens that the developers of [IBDMDB](https://www.ibdmdb.org/) have provided the pre-processed data to researchers so it would have been strange to find many low-quality reads. 

**Question 3: It's reasonable that the processed data was provided by IBDMDB since their goal is to provide standardized data that can be used for many projects. However, why is it important that _normally_ all raw data be uploaded to public repositories?**

Let's enter the kneaddata output folder now and see how large the output files are:

```
cd kneaddata_out
ls -lh *fastq
```

In this example case the output files are quite small (~5 megabytes), but for typical datasets these files can be extremely large so you may want to delete the intermediate FASTQs if you're running your own data through this pipeline. Also, note that the output FASTQs are not gzipped - you could gzip these files if you planed on saving them long-term. In the meantime you can move the intermediate FASTQs to their own directory with these commands:

```
mkdir intermediate_fastqs
mv *_contam*fastq intermediate_fastqs
mv *unmatched*fastq intermediate_fastqs
cd ..
```
### Processing forward and reverse reads

We have one last step before we are ready to begin generating taxonomic profiles! At this point we have reached a fork in the road where we have two options. The first option is to stitch our forward and reverse reads together while the second is to simple concatenate (combined) the forward and reveres reads into a single sequencing file.

#### Stitching reads

The process of stitching involves taking the forward read and combining it with the reverse complement of the reverse read based on where the reads would overlap in the middle.

This process can result in the majority of your reads being thrown out during the stitching process as many metagenomic reads will have little to no overlap. This is because in many cases during library preparation the DNA fragments are too large for the forward and reverse reads to completely traverse them. Remember that the forward and reverse reads only cover 100-150bp depending on platform and DNA fragments from many metagenomic library preps average in size between 300-500bp. This means there is not sufficient overlap to bring the forward and reverse reads together. As such I would suggest skipping this step to avoid losing a majority of reads.

If you would like to read more about how you might go about the stitching process using the program [PEAR](https://academic.oup.com/bioinformatics/article/30/5/614/247231) check out this page [here](https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-standard-operating-procedure-v2#1-first-steps).

#### Concatenating reads
Instead of dealing with the issue of stitching we will instead concatenate (combined) the forward and reverse reads into a single file. This is why its important during our filtering steps to only keep reads that both pairs passed all filtering parameters. If not we might end up with unaccounted for bias due to more forward reads being included in our concatenated reads than reverse. 

We can do this by running the command:


```
ln -s ~/CourseData/MIC_data/Module4/Data/helper_scripts/ .
perl helper_scripts/concat_paired_end.pl -p 4 --no_R_match -o cat_reads kneaddata_out/*_paired_*.fastq
```

Note the first command is simply linking the folder that contains the script for running this command to the current directory we are working in. 


***Always double check the output to make sure the correct files are being concatenated together***

## Generating Taxonomic Profiles with Kraken2 + Bracken

Now that we have processed all of our reads the next step is to use [Kraken2](https://ccb.jhu.edu/software/kraken2/) to classify them with taxonomic labels. Kraken2 is one of many popular tools for classifying metagenomic reads into taxonomic profiles. In this case Kraken2 uses a kmer approach to assign reads to specific taxonomic lineages (look back at the lecture slides to remind yourself how kmers are used to classify reads). 

### Databases
Unlike metagenomic assembly we will be classify our sequenced reads based on how similar they are to reads we have already classified in a reference database. Kraken2 comes with a default database that is 32Gb, much too large to be run on our instances that only have 16Gb of RAM/memory. Therefore we will have to use the a secondary database supplied by the Kraken2 develops called [minikraken2](http://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads). This database is only 8Gb in size which will allow you to run it on workstations with a lower amount of RAM at the cost of reduced performance. Note that by default kraken uses NCBI taxonomy, however, there are a number of other databases that users have made using different styles of taxonomy such as the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/)

We would not normally recommend the use of these smaller databases such as minikraken2 unless the users cannot get access to a workstation with a large amount of RAM. This is because these smaller databases can suffer from lower levels of classification especially if you increase the confidence threshold used by kraken2. 


### Running Kraken2 + Bracken
We are going to use our knowledge about parallel from above to process both of our samples through kraken2 one after another!

We will first make directories to hold our results
```
mkdir kraken2_outraw
mkdir kraken2_kreport
```

Next we will run the Kraken2 command using parallel.

```
parallel --dry-run -j 1 'kraken2 --db  ~/CourseData/MIC_data/Module4/Data/kraken_db/minikraken2_v2_8GB_201904_UPDATE/ --threads 4 --output kraken2_outraw/{/.}.kraken.txt --report kraken2_kreport/{/.}.kreport --use-names {}' ::: cat_reads/*.fastq
```
* ```--db``` - This option points to the location of the database you want Kraken2 to use for classification. We have already set this AWS instance with the minikraken database.
* ```--threads 4```-This option indicates that we want to use four processors for each Kraken2 job
* ```--output``` - This argument points to the file name we want to output our classifications to.
* ```--report``` - This argument points to the file name we want to output a more detailed classifcation report for each sample. These more in-depth classification reports are required to correct abundance estimates using bracken. 
* ```--use-names``` - This indicates we want the classifier to output the taxon names rather than their NCBI IDs. 

**Question 4: What does the text ```{/.}.kraken.txt``` and ```{/.}.kreport``` do in the above command? (Hint look at our above discussion of parallel).**

Check to make sure the command you are running is the one you expected using the --dry-run flag. Once you are sure remove the flag and run the command.

Some summary information will be output to the console when you run this command including the number of sequences that were classified versus unclassified.

#### Correcting raw kraken2 output with bracken
The raw output files created by Kraken2 should be corrected using Bracken to get more accurate estimations on the abundance of different taxa at various taxonomic levels. While Kraken2 classifies each of our reads to their best locations within a taxonomic tree. In some cases this will result in a read being placed higher up in the tree due to two species or sometimes even genera sharing the exact same sequence (kmer) to which that read matched. Bracken can be used to solve this issue by taking these reads and mapping them back to the genomes they might belong to. By doing this it generates new abundance estimates that take into account the problem laid out above.

To apply Bracken to our samples to estimate the species level abundances we can run the following commands:

First we will make a directory to output the results too:

```
mkdir bracken_out
```
Now we will run the bracken command:

```
source activate CBW_2021_MARKER
parallel --dry-run -j 2 'bracken -d ~/CourseData/MIC_data/Module4/Data/kraken_db/minikraken2_v2_8GB_201904_UPDATE/ -i {} -o bracken_out/{/.}.species.bracken -r 100 -l S -t 1' ::: kraken2_kreport/*.kreport
source activate CBW_2021_KMER
```
* ```-i ``` -The input kreport file from kraken2
* ```-r``` - The length of the reads used for mapping
* ```-l``` - The level at which to estimate taxon abundance (S indicates species)
* ```-t``` - The number of reads required for a taxon to be used in the abundance estimate

Note generally you would increase the -t argument to be somewhere between 0.1-1% of the total average read depth however since we are working with files that have already been subset we will not set this option. The reason for this is that kmer based strategies such as kraken2 are notorious for making a large number of low abundance false positives. Check-out this [paper](https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5) for more information.

The final step in generating our taxonomic profiles is collecting all the results from each sample into a single file. This can be done using a helper script that comes installed with bracken!

```
mkdir bracken_out_merged
python helper_scripts/combine_bracken_outputs.py --files bracken_out/*.species.bracken -o bracken_out_merged/merged_output.species.bracken
```

Lets take a look at our resulting file with the ```less``` command. You will notice that the file has 7 different columns corresponding to different information about a single taxon.

| Column Name      | Data Type     |
| :------------- | :----------: |
| name |  Taxonomic ID  |
| taxonomy_id | NCBI taxonomic ID |
| taxonomy lvl | Character representing the level of taxonomy |
| SAMPLE.species.bracken_num | Number of reads that aligned to this classifcation for this SAMPLE|
| SAMPLE.species.bracken_frac | Proportion of reads that aligned to this classification for this SAMPLE |
| SAMPLE2.species.bracken_num | Number of reads that aligned to this classifcation for this SAMPLE2 |
| SAMPLE2.species.bracken_frac | Proportion of reads that aligned to this classification for this SAMPLE2 |


***Question 5: How would you go about figuring out the total number of species we found within our two different samples?***

## Determining functional profiles using MMseqs2

Now that we have an idea of the taxonomic profiles of our samples we can also look at their functional potential. To do this we will being using [MMseqs2](https://github.com/soedinglab/mmseqs2) along with a few scripts developed by our lab to connect the functional annotations to the taxonomic annotations.

MMseqs2 works by taking sequenced reads translating them into protein and then mapping them against a protein database. In this case we will be using [UniRef90](https://www.uniprot.org/help/uniref). A large protein database clustered at 90% identity. Our lab has developed a set of scripts to run this tool and use the information to assign both a taxonomic ID and functional ID to each read within our samples.

### Job file generation

One way of running mutliple commands in terminal one after another is to create a job file. This job file contains all of the commands of interest along with all the needed parameters to run those commands. We can use the following command to create a job file for each sample that we would like to process with MMseqs2.

```
ln -s ~/CourseData/MIC_data/Module4/Data/Functional_Helper_Scripts/ .
/usr/bin/perl Functional_Helper_Scripts/run_TaxonomyFunctionSearchMegan.pl --mapmethod mmseqs \
--db MMSeqs2_db/mmseqsUniref90DB -p 4 -o mmseqs_U90_out cat_reads/*fastq
```

This should create two separate job files one for each sample. Well you don't need to look in these files its always a good idea to understand what your bioinformatic programs are doing! Lets take take a look at one of the job files using the ```less``` command. You should see three lines each containing a separate command.

```
mmseqs createdb cat_reads/CSM79HR8.fastq mmseqs_U90_out/mmseqs-CSM79HR8queryDB
```

This command creates an mmseqs database from the the input fastq file. The creation of this database is necessary for mmseqs as it vastly increases the speed at which translated DNA sequences can be mapped against a protein database.

```
mmseqs search mmseqs_U90_out/mmseqs-CSM79HR8queryDB MMSeqs2_db/mmseqsUniref90DB mmseqs_U90_out/mmseqs-CSM79HR8resultDB tmp --db-load-mode 3 --threads 4 --max-seqs 25 -s 1 -a -e 1e-5 > /dev/null 2>&1
```

This command is the real meat of the job file and runs the freshly created sample database against the provided UniRef90 protien database. There are a number of parameters in this command.
* ```--db-load-mode 3``` - This parameter tells mmseqs how to deal with loading the database into memory. For more information you can check out this [page]. However, setting this parameter to 3 helps when running mmseqs on a cluster environment. 
* ```--threads``` - The number of processors we want mmseqs to use during the search
* ```--max-seqs 25``` - This indicates that we want mmseqs to output at maximum 25 hits for each sequence
* ```-s 1``` - This indicates the sensitivity that we want mmseqs to run at. Increasing this number will lower the speed at which mmseqs runs but will increase its sensitivity. For well explored environments a setting of 1 should suffice.
* ```-a``` - Adds indicates that we want our results to output backtraces for each sequence match. These are needed to convert the resulting mmseqs file into a usable file format. 
* ```-e 1e-5``` - This indicates that we only want to keep matches that are below an E-value of 1e-5 (E-values are a measure of how well two sequences match one another).
* ```> /dev/null 2>&1``` This final part of the command allows us to run the command without having too much text printed to our screen.

The final command allows us to convert the resulting file from the MMSeqs2 format into one that is more usable.

```
mmseqs convertalis mmseqs_U90_out/mmseqs-CSM79HR8queryDB MMSeqs2_db/mmseqsUniref90DB mmseqs_U90_out/mmseqs-CSM79HR8resultDB mmseqs_U90_out/mmseqs-CSM79HR8-s1.m8 --db-load-mode 2 > /dev/null 2>&1
```

This command is similar and takes as input the query database we made from our first command, the UniRef90 database we searched against and the resulting file from our search command. It will output the file ```mmseqs_U90_out/mmseqs-CSM79HR8-s1.m8```. 

**Unfortunately the AWS instances we are using do not have enough memory to run the above three commands in a timely manner (they over an hour on this AWS configuration). As such I already provided the expected output from the above commands. You can copy them to your working directory using the following commands:**

```
cp -r ~/CourseData/MIC_data/Module4/Intermediate_Results/Module4/mmseqs_U90_out .
cp -r ~/CourseData/MIC_data/Module4/Intermediate_Results/Module4/mmseqs_m8_files .
```

All of these commands made a lot of files that we are actually not interested in! We will go ahead and remove all of the files except for those that end in .m8. This will help us save significantly on the amount of hard drive space that is being taken up by these files.

```
mkdir mmseqs_m8_files
mv mmseqs_U90_out/*.m8 mmseqs_m8_files/
rm -r mmseqs_U90_out
```

**The above 3 commands have already been run for you.**

Lets take a quick look at one of the files we just moved into the directory ```mmseqs_m8_files``` using the less command.

```
less mmseqs_m8_files/mmseqs-CSM79HR8-s1.m8
```

We you will see is a file in BLAST tabular format.

| Column Number       | Data Type     |
| :------------- | :----------: |
| 0 |  query sequence ID  |
| 1 | Subject (database) sequence ID |
| 2 | 	Percent Identity |
| 3 | Alignment Length |
| 4 | Number of gaps |
| 5 | 	Number of mismatches |
| 6 | Start on the query sequence |
| 7 | End on the query sequence |
| 8 | Start on the database sequence |
| 9 | 	End on the database sequence |
| 10 | 	E value - the expectation that this alignment is random given the length of the sequence and length of the database |
| 11 | bit score - the score of the alignment itself |


**Question 6: How many protein sequences did the sequence ```HKWJVBCXY170606:2:2116:8029:10262/1``` align with in the sample CSM79HR8? What alignment/alignments have the lowest E-value/highest bitscore?**


The next step we need to take is to get the name of the protein sequence that had the best alignmnet for each sequence read in our samples. We can achieve this by running the command:


```
mkdir mmseqs_U90_out_tophit
python Functional_Helper_Scripts/pick_uniref_top_hit.py --unirefm8Dir mmseqs_m8_files --output_path mmseqs_U90_out_tophit
```
Now that we have the best protein sequence that matches best with each of the sequences in our samples we can begin creating our final data table containing the stratified abundance of each function in our samples.

First we need to generate an input file that contains information about where all of the results we have generated are located. We can generate this table using the following commands:

* The 1st step is to generate the sample tags; we can use either the concatenated fastq files or the kraken output for extracting the sample tag IDs

```ls -1 kraken2_outraw/*kraken.txt |sed -e 's/.*\///g'|sed -e 's/\.kraken\.txt//g' > sample-tags.txt```

* Next we need the paths to the uniref funtional files (parsed to get only the tophit) and the filetype (refseq, uniref or COG)

```ls -1 mmseqs_U90_out_tophit/* | sed -e 's/$/\tuniref/g' > uniref-hits-list.txt```

* Similarly, we would need the paths to the kraken2 classification files and the corresponding filetypes

```ls -1 kraken2_outraw/*.kraken.txt | sed -e 's/$/\tkraken2/g' > kraken2-results-list.txt```

* Finally, we need the list of the m8 output files for each sample

```ls -1 mmseqs_m8_files/* > m8-list.txt```

* The last step is to paste put all of these files together

```paste sample-tags.txt kraken2-results-list.txt uniref-hits-list.txt m8-list.txt > multi-sample-outfiles-w-m8.txt```

Now that we have this master file we can pass this information into the 

```
ln -s ~/CourseData/MIC_data/Module4/Data/MMSeqs2_db/*.pbz2 .
python Functional_Helper_Scripts/parse_TaxonomyFunction.py --multisample multi-sample-outfiles-w-m8.txt --outputf Workshop_strat_matrix_RPKM.txt --stratified Y --map2EC Y
```

The first link command will grab all the databases that contain information about the length of each gene in our UniRef protein database. This will be important to normalize the abundance of each functional classification. 

The second command will generate a final stratified data matrix that shows the abundance of each EC number stratified by the different taxonomic classifications within the sample. This script also normalizes the abundances of each of these ECS into reads per kilobase per million (RPKM). This abundance metric is the number of reads that mapped to that EC number per million reads mapped within the sample divided by the gene length of that EC. Its important that the abundances are normalized by gene length or there would be an unfair bias toward longer genes due to the higher chance of them being sequenced. We can also run the same command as above without the ```--stratified Y``` option to generate functional profiles that are not broken down by the contributing taxa. 

**Wait a minute my command ended with the output ```Killed```!**

Do not panic this was expected as the AWS instance we are running on does not have enough memory to load the database that contains the lengths of each gene in our UniRef90 database. However I have already generated the expected output from this command and we can copy it over to our working directory!

```
cp ~/CourseData/MIC_data/Module4/Intermediate_Results/Module4/Workshop_strat_matrix_RPKM.txt .
```

Lets take a look at this file with the less command:

```
less Workshop_strat_matrix_RPKM.txt
```
*** Question 7: What is RPKM contributed to the sample CSM79HR8 for the EC 2.1.2.9 contributed by Bacteroides vulgatus? ***


### Adding descriptions and calculating pathway abundances using PICRUSt2


Using scripts from PICRUSt2 we can add descriptions (i.e the names for each EC number) or generate stratified pathway abundances. Before we can do this we need to make some minor formatting adjustments.

```
sed -e "s/|/\t/g" Workshop_strat_matrix_RPKM.txt > Workshop_strat_matrix_RPKM_fixed.txt
sed -e "s/taxonomy/sequence/g" Workshop_strat_matrix_RPKM_fixed.txt > Workshop_strat_matrix_RPKM_fixed_final.txt
```

Now we will activate the picrust2 conda environment.

```
source activate picrust2
```

We can add descriptions using the following command: 

```
add_descriptions.py -i Workshop_strat_matrix_RPKM_fixed_final.txt -m EC -o Workshop_strat_matrix_RPKM_fixed_des.txt
```

**Question 8: What is the name of the enzyme with the EC number 6.1.1.4?**

We can generate stratified pathway abundances using the following command:

```
pathway_pipeline.py -i Workshop_strat_matrix_RPKM_fixed_final.txt -p 4 -o pathways_out --wide_table
```

The results of this command will be in the pathways_out folder. Just like in the picrust2 tutorial! 

**Question 9: How many total pathways were identified?**

**Question 10: how many taxa contribute to the pathway ```PANTO-PWY```?**

## Next Steps to try if we have time!

Now that we have generated data tables that include the taxonomic breakdowns of samples using Kraken2+Bracken and the stratified functional contributions using MMseqs2 and Kraken2 we can now use various tools to analysis the data. Tools that you might be interested in using to view the stratified functional contributions include [BURRITO](http://borenstein-lab.github.io/burrito/) or making use of the pathway tables using [Metacyc's Omics Dashboard](https://www.metacyc.org/dashboard/dashboard-intro.shtml).

If your interesting in playing around with these tools and would like more than two samples you can download processed taxonomic and functional tables for 18 additional samples using the following commands.

Bracken genus level profiles:
```
wget https://www.dropbox.com/s/dqhgnydrnyrugy7/merged_output.genus.bracken?dl=0
```

Stratified functional contributions:
```
wget https://www.dropbox.com/s/oe898c7uqn0lrpy/FULLDATA_strat_matrix_RPKM_fixed.txt?dl=0
```

The metadata from these samples can be moved into your working directory using these commands.


```
cp ~/CourseData/MIC_data/Module4/Data/mgs_metadata.txt .
```

If you have time I have prepared an additional R script that you can run to create a stacked bar chart out of the taxonomic profiles we just downloaded.
You can download and run this script using the following commands:

```
wget https://www.dropbox.com/s/ccc3a8tiqi9xcnj/Bracken_Stacked_Barcharts.R?dl=0
mv Bracken_Stacked_Barcharts.R\?dl\=0 Bracken_Stacked_Barcharts.R
Rscript Bracken_Stacked_Barcharts.R
```

You can now open your instance in your browser and open up the PDF ```Taxa_barplots.pdf``` that was generated from running this script. If you have time I would also suggest inspecting the Rscript for an example of how you might load bracken data into R and create visualizations such as stacked barcharts. 
