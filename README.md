# STAR-docker
[![Docker Image CI](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-image.yml/badge.svg)](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-image.yml)[![Docker](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-publish.yml)

A minimal container holding [STAR 2.7.10b](https://github.com/alexdobin/STAR).
Built from the [GitHub Release](https://github.com/alexdobin/STAR/releases/tag/2.7.10b) using wget, base image off of ubuntu

[STAR is now Splice aware](https://github.com/alexdobin/STAR#readme)

[Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2022](https://www.ncbi.nlm.nih.gov/pubmed/23104886)

[STAR User Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

### Requirements

This was built on a Mac Book Pro.
Docker is installed -- need root privileges for installation


### Getting the needed bits and parts

Actually grabbing a couple of small test data sets from our [Dry Bench Skills Class](https://doi.org/10.5281/zenodo.7025773)

We set up a test directory after cloning this repository for testing.

To test out this installation and docker build which will be stitched together into our workflow, we use test data From a previous workflow for testing our Proteogenomics Long Read RNA-Sequencing Workflow [Test Data](https://doi.org/10.5281/zenodo.5234651)

From here we will get our genome and annotation small files for testing

You can test the container by doing the following:

First, make a test directory (just to keep things clean)
```bash
mkdir test
cd test
```

Then download two test files
```bash
wget https://zenodo.org/record/7025773/files/test.20k_reads_1.fastq.gz
wget https://zenodo.org/record/7025773/files/test.20k_reads_2.fastq.gz
```

Now download a genome with annotations (actually just chromosome 22 -- to keep it small.

```bash
wget https://zenodo.org/record/5234651/files/GRCh38.primary_assembly.genome.chr22.fa
wget https://zenodo.org/record/5234651/files/gencode.v35.annotation.chr22.gtf
```

Next, referring to the provided STAR manual we first index the genome

```bash
docker build -t star .
```

To test this tool from the command line 

Set up an environment variable capturing your current command line:
```bash
PWD=$(pwd)
```

Then mount and use your current directory and call the tool now encapsulated within the environment.


## Test your docker image

First, let's test to see the docker image was successful

```bash
docker run -it -v $PWD:$PWD -w $PWD star STAR --version
```

## Now lets build an index

Referring to the manual, we see we can generate an index with the *genomeGenerate* command

```bash
docker run -it -v $PWD:$PWD -w $PWD star STAR \
     --runMode genomeGenerate \
     --runThreadN=7 \
     --genomeFastaFiles GRCh38.primary_assembly.genome.chr22.fa \
     --genomeDir . \
     --sjdbGTFfile gencode.v35.annotation.chr22.gtf \
     --sjdbOverhang 50 \
     --outFileNamePrefix chr22
```

This was run on a Macbook Pro (see configuration of the machine) taking about 1 minute with chr22 only file.

<p>
<br/><br/>
<img src="https://github.com/adeslatt/STAR-docker/blob/main/assets/MacBookProDetails2022Dec05.png" width=500">
<br/><br/>
</p>

If successful, the *`Log.out`* file resulting from running in *`runMode`* to get the indices looks like this:

```bash
STAR version=2.7.10b
STAR compilation time,server,dir=2022-11-01T09:53:26-04:00 :/home/dobin/data/STAR/STARcode/STAR.master/source
STAR git: On branch master ; commit c6f8efc2c7043ef83bf8b0d9bed36bbb6b9b1133 ; diff files: CHANGES.md 
##### Command Line:
STAR --runMode genomeGenerate --genomeFastaFiles GRCh38.primary_assembly.genome.chr22.fa --genomeDir . --sjdbGTFfile gencode.v35.annotation.chr22.gtf --sjdbOverhang 50 --outFileNamePrefix chr22
##### Initial USER parameters from Command Line:
outFileNamePrefix                 chr22
###### All USER parameters from Command Line:
runMode                       genomeGenerate        ~RE-DEFINED
genomeFastaFiles              GRCh38.primary_assembly.genome.chr22.fa        ~RE-DEFINED
genomeDir                     .     ~RE-DEFINED
sjdbGTFfile                   gencode.v35.annotation.chr22.gtf     ~RE-DEFINED
sjdbOverhang                  50     ~RE-DEFINED
outFileNamePrefix             chr22     ~RE-DEFINED
##### Finished reading parameters from all sources

##### Final user re-defined parameters-----------------:
runMode                           genomeGenerate   
genomeDir                         .
genomeFastaFiles                  GRCh38.primary_assembly.genome.chr22.fa   
outFileNamePrefix                 chr22
sjdbGTFfile                       gencode.v35.annotation.chr22.gtf
sjdbOverhang                      50

-------------------------------
##### Final effective command line:
STAR   --runMode genomeGenerate      --genomeDir .   --genomeFastaFiles GRCh38.primary_assembly.genome.chr22.fa      --outFileNamePrefix chr22   --sjdbGTFfile gencode.v35.annotation.chr22.gtf   --sjdbOverhang 50
----------------------------------------

Number of fastq files for each mate = 1
ParametersSolo: --soloCellFilterType CellRanger2.2 filtering parameters:  3000 0.99 10
Finished loading and checking parameters
--genomeDir directory exists and will be overwritten: ./
Nov 10 23:47:12 ... starting to generate Genome files
GRCh38.primary_assembly.genome.chr22.fa : chr # 0  "chr22" chrStart: 0
Chromosome sequence lengths: 
chr22   50818468
Genome sequence total length = 50818468
Genome size with padding = 50855936
Nov 10 23:47:13 ..... processing annotations GTF
Processing pGe.sjdbGTFfile=gencode.v35.annotation.chr22.gtf, found:
                4964 transcripts
                30205 exons (non-collapsed)
                8165 collapsed junctions
Total junctions: 8165
Nov 10 23:47:14 ..... finished GTF processing

!!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=50818468, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 11
Estimated genome size with padding and SJs: total=genome+SJ=151855936 = 50855936 + 101000000
GstrandBit=32
Number of SA indices: 78319554
Nov 10 23:47:14 ... starting to sort Suffix Array. This may take a long time...
Number of chunks: 1;   chunks size limit: 18538972752 bytes
Nov 10 23:47:14 ... sorting Suffix Array chunks and saving them to disk...
Writing 626556432 bytes into .//SA_0 ; empty space on disk = 225814350135296 bytes ... done
Nov 10 23:47:52 ... loading chunks from disk, packing SA...
Nov 10 23:47:56 ... finished generating suffix array
Nov 10 23:47:56 ... generating Suffix Array index
Nov 10 23:48:13 ... completed Suffix Array index
Nov 10 23:48:13   Finished preparing junctions
Nov 10 23:48:13 ..... inserting junctions into the genome indices
Nov 10 23:48:20   Finished SA search: number of new junctions=8162, old junctions=0
Nov 10 23:48:21   Finished sorting SA indicesL nInd=1632400
Genome size with junctions=51680298  50855936   824362
GstrandBit1=32   GstrandBit=32
Nov 10 23:48:22   Finished inserting junction indices
Nov 10 23:48:28   Finished SAi
Nov 10 23:48:28 ..... finished inserting junctions into genome
Nov 10 23:48:28 ... writing Genome to disk ...
Writing 51680298 bytes into .//Genome ; empty space on disk = 225262722613248 bytes ... done
SA size in bytes: 329801814
Nov 10 23:48:29 ... writing Suffix Array to disk ...
Writing 329801814 bytes into .//SA ; empty space on disk = 225249490632704 bytes ... done
Nov 10 23:48:32 ... writing SAindex to disk
Writing 8 bytes into .//SAindex ; empty space on disk = 225437401743360 bytes ... done
Writing 120 bytes into .//SAindex ; empty space on disk = 225437404889088 bytes ... done
Writing 1565873491 bytes into .//SAindex ; empty space on disk = 225437424812032 bytes ... done
Nov 10 23:48:41 ..... finished successfully
DONE: Genome generation, EXITING
```
## (Optional) Build your genome index in a Terminal Window within a JupyterLab notebook within [CAVATICA](cavatica.sbgenomics.com) 

So wait a minute... I don't have enough RAM on my 16 GB RAM MacBook Pro, what am I thinking???

Ok, so what is the work around.   Can't run docker within a terminal window within a JupyterLab notebook on CAVATICA.  Security problem apparently.  FedRAMP issue.  So what to do.

### Anaconda install STAR

Conda is your friend.  Conda is on CAVATICA and may be found within the image of the environment, and because we are in the cloud you can spin up as a large a machine deemed necessary for your genome.  Remember not all genomes are the same.  Trisomy 21 for example.  What is that reference genome.  I'll get back to you on that one.  

How do you do it?

#### [Anaconda search packages star](anaconda.org/bioconda/star)

Nearly 399,918 downloads, I believe you can trust this install.  Bioconda, look it up, a trust-worthy packager of packages and in many environments, also here in conda.  The bioconda channel - again look it up, you may learn something -- keeps it up todate.   Inspecting you see that it is up-to-date with the *`v2.7.10b`* that we want to use.  Both on the Mac (osx-64), as well as the linux-64.  

Navigate to a terminal window within the `*Jupyterlab`* Launcher and install **`star`**

```bash
conda install -c bioconda star
```

This will install **`star`**

```bash
STAR -h
```

Will spit out a lot of help information.

#### Get reference genome and annotation

For this exercise, we need gencode annotation v33.   Googling it we find where to grab the appropriate genome input fasta file and annotation files

For clean directory set up (we will later move this to the project files on the CAVATICA platform, but for now we will work in a subdirectory to get the work done.

```bash
mkdir human_reference
cd human_reference
```

Now, lets get our dependencies: a genome fasta file and an annotation file.

First the annotation file:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
gunzip gencode.v33.annotation.gtf.gz
```

**`STAR`** can deal with gzipped input fastq files, but not annotation and genome files

Now the genome reference file.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

Now we are ready to build our genome indices.  Once we have them, we will move them to the project-files location for all of us to take advantage of.


```bash
STAR \
--runMode genomeGenerate \
--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v33.annotation.gtf \
--sjdbOverhang 50  \
--runThreadN=14    \
--outFileNamePrefix gencode.v33.
```

Note that I spun up a 16 vCPU with 32 GB RAM.  Note that most often you are getting about 2GB RAM per vCPU.  Just good to know.

## Test with test reads

First, let's test to see if the java install in the docker image was successful

Using some test reads from the same repository where the above chr22 limited genome was found.
For the curious, it came from course I taught with the wonderful Christina Chatzipantsiou @cgpu.
[Dry Bench Skills for the Researcher](https://doi.org/10.5281/zenodo.7025773).   This was done with [Lifebit.ai](https://www.lifebit.ai/) while we were at [The Jackson Laboratory](https://www.jax.org)
                                                         
```bash
wget https://zenodo.org/record/7025773/files/test.20k_reads_1.fastq.gz
wget https://zenodo.org/record/7025773/files/test.20k_reads_2.fastq.gz
```
Referring to the [STAR User Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) we see how we can specify
* How to specify our input fastq files (note that ours are compressed)
* How to specify number of threads (on the mac I'll choose 6) to enhance parallel processing
* How to input our reads (we could use a manifest file)
* How to specify we want coordinate sorted BAM file (I had to increase on the mac the number of files that can be created **`ulimit -n 5000`** - the default on the mac was **`256`**
* How to specify annotation files (we will use gencode.v35.annotation.chr22.gtf as we did with our genome)
* How to specify plicing junction list (sjdbList.fromGTF.out.tab output in the genome index creation step above)
* How to provide a outFileNamePrefix
* How to output header details for the bam (outSAMheader @HD, outSAMheaderPG @PG and outputSAMheaderCommentFile @CO)
* and how to ensure BAM sorting uses the number of threads available (outBAMsortingThreadN=6)
                                                         
Now we are ready to go                                              
       
```bash
docker run -it -v $PWD:$PWD -w $PWD star STAR \
     --runThreadN=6 \
     --genomeDir $PWD \
     --readFilesIn /Users/deslattesmaysa2/projects/STAR-docker/test/test.20k_reads_1.fastq.gz \
                   /Users/deslattesmaysa2/projects/STAR-docker/test/test.20k_reads_2.fastq.gz \
     --sjdbGTFfile /Users/deslattesmaysa2/projects/STAR-docker/test/gencode.v35.annotation.chr22.gtf \
     --sjdbFileChrStartEnd /Users/deslattesmaysa2/projects/STAR-docker/test/sjdbList.fromGTF.out.tab \
     --outFileNamePrefix $PWD/test.chr22. \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand gunzip -c \
     --outSAMheaderHD @HD \
     --outSAMheaderPG @PG \
     --outSAMheaderCommentFile @CO \
     --outBAMsortingThreadN=6
```

Runs very quickly on the Mac Book Pro
Successful execution looks like this:

```bash
STAR --runThreadN=6 --genomeDir /Users/deslattesmaysa2/projects/STAR-docker/test --readFilesIn /Users/deslattesmaysa2/projects/STAR-docker/test/test.20k_reads_1.fastq.gz /Users/deslattesmaysa2/projects/STAR-docker/test/test.20k_reads_2.fastq.gz --sjdbGTFfile /Users/deslattesmaysa2/projects/STAR-docker/test/gencode.v35.annotation.chr22.gtf --sjdbFileChrStartEnd /Users/deslattesmaysa2/projects/STAR-docker/test/sjdbList.fromGTF.out.tab --outFileNamePrefix /Users/deslattesmaysa2/projects/STAR-docker/test/test.chr22.withHeader. --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --outSAMheaderHD @HD --outSAMheaderPG @PG --outSAMheaderCommentFile @CO --outBAMsortingThreadN=6
	STAR version: 2.7.10b   compiled: 2022-11-01T09:53:26-04:00 :/home/dobin/data/STAR/STARcode/STAR.master/source
Dec 11 22:00:31 ..... started STAR run
Dec 11 22:00:31 ..... loading genome
Dec 11 22:00:40 ..... processing annotations GTF
Dec 11 22:00:41 ..... inserting junctions into the genome indices
Dec 11 22:00:45 ..... started mapping
Dec 11 22:02:23 ..... finished mapping
Dec 11 22:02:23 ..... started sorting BAM
Dec 11 22:02:24 ..... finished successfully
```

And with the parameters above produces 5 output files:

```bash
-rwxrwxr-x  1 deslattesmaysa2  NIH\Domain Users        9526 Dec 11 16:37 test.chr22.Log.out
-rw-r--r--  1 deslattesmaysa2  NIH\Domain Users        1981 Dec 11 16:37 test.chr22.Log.final.out
-rwxrwxr-x  1 deslattesmaysa2  NIH\Domain Users         364 Dec 11 16:37 test.chr22.Log.progress.out
-rwxrwxr-x  1 deslattesmaysa2  NIH\Domain Users      125902 Dec 11 16:37 test.chr22.Aligned.sortedByCoord.out.bam
-rw-r--r--  1 deslattesmaysa2  NIH\Domain Users         446 Dec 11 16:37 test.chr22.SJ.out.tab
```

## (Optional) Deposit your container in the your [CAVATICA](cavatica.sbgenomics.com)  Docker Registry

If you are working with Kids First or INCLUDE data, and you are registered with either the [Kids First DRC](https://kidsfirstdrc.org/) or the [INCLUDE Data Hub](https://includedcc.org/), you have access to a platform as a service, CAVATICA by Seven Bridges.

If you do, you Docker Image Repository Specific to you is at the location of pgc-images.sbgenomics.com/[YOUR CAVATICA USERNAME]

There are three steps to building and using a container

1. build
2. tag
3. push

We have built and we did tag this (I built the above on my desktop which is a mac, but it could have been built within a google shell.  Examples of how to do this are found in the [Elements of Style in Workflow Creation and Maintenance, Day 3](https://github.com/NIH-NICHD/Kids-First-Elements-of-Style-Workflow-Creation-Maintenance/blob/main/classes/Building-A-Nextflow-Script/README.md#preamble-to-building-workflows-using-containers)

### Tag

To tag the image just built, you need the image id, to get that simply use the command **`docker images`**.

```bash
docker images
```

This returns:

```bash
REPOSITORY                                           TAG       IMAGE ID       CREATED        SIZE
samtools                                             latest    e7945dfbaa32   14 hours ago   576MB
rmats                                                latest    24f836d5d5a2   6 days ago     1.53GB
star                                                 latest    264b10d079ae   4 weeks ago    125MB
trimmomatic                                          latest    633fec609229   4 weeks ago    790MB
```

Now we re-tag it for pushing to our own personal [CAVATICA](cavatica.sbgenomics.com) docker container registry.

```bash
docker tag 264b10d079ae pgc-images.sbgenomics.com/[YOUR CAVATICA USERNAME]/star:v2.7.10b
```

### Docker registry login

There is actually another step required before you can do your push to the registry.  You need to authenticate.

Navigate to the CAVATICA Developers Tab.

Select Authentication Token, if you have not done so, generate that token.

Then Copy the token and paste it in the proper location in the command below (after -p for password)

```bash
docker login pgc-images.sbgenomics.com -u [YOUR CAVATICA USERNAME] -p [YOUR AUTHENTICATION TOKEN]
```

### Push

Now that we have

* tagged :whitecheck our docker image and

* we have authenticated :whitecheck 

```bash
docker push pgc-images.sbgenomics.com/[YOUR CAVATICA USERNAME]/star:v2.7.10b
```

You know things are going correctly when you see something to the effect of:


```bash
The push refers to repository [pgc-images.sbgenomics.com/[YOUR CAVATICA USERNAME]/star:v2.7.10b]
ea25457229f1: Pushed 
17280cc0fa6b: Pushed 
ab2731ec3f53: Pushed 
6fa1f4185aa2: Pushed 
ad6562704f37: Pushed 
latest: digest: sha256:3b1976baa7c4aaa2afd25098c41c754e2579060b6c1da32282c45ac8a10293a9 size: 1373
```
