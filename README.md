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

First, let's test to see if the java install in the docker image was successful

```bash
docker run -it -v $PWD:$PWD -w $PWD star STAR --version
```

## Now lets build an index

Referring to the manual, we see we can generate an index with the *genomeGenerate* command

```bash
docker run -it -v $PWD:$PWD -w $PWD star STAR \
     --runMode genomeGenerate \
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



