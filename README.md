# STAR-docker
[![Docker Image CI](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-image.yml/badge.svg)](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-image.yml)[![Docker](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/adeslatt/STAR-docker/actions/workflows/docker-publish.yml)

A minimal container holding [STAR 2.7.10b](https://github.com/alexdobin/STAR).
Built from the [GitHub Release](https://github.com/alexdobin/STAR/releases/tag/2.7.10b) using wget, base image off of ubuntu

[STAR is now Splice aware](https://github.com/alexdobin/STAR#readme)

[Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2022](https://www.ncbi.nlm.nih.gov/pubmed/23104886)

[STAR User Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

## Using Test Data

To test out this installation and docker build which will be stitched together into our workflow, we use test data From a previous workflow for testing our Proteogenomics Long Read RNA-Sequencing Workflow [Test Data](https://doi.org/10.5281/zenodo.5234651)

### Getting the needed bits and parts

Actually grabbing a couple of small test data sets from our [Dry Bench Skills Class](https://doi.org/10.5281/zenodo.7025773)

We set up a test directory after cloning this repository for testing.

#### Requirements

This was built on a Mac Book Pro.
Docker is installed -- need root privileges for installation

#### Small Fastq data sets

Then download two test files
```bash
wget https://zenodo.org/record/7025773/files/test.20k_reads_1.fastq.gz
wget https://zenodo.org/record/7025773/files/test.20k_reads_2.fastq.gz
```

Next, referring to the provided STAR manual we first index the genome

```bash
docker build -t trimmomatic .
```

To test this tool from the command line 

Set up an environment variable capturing your current command line:
```bash
PWD=$(pwd)
```

Then mount and use your current directory and call the tool now encapsulated within the environment.

Any trimmomatic is a java application and running it requires the jar file which is unzipped in this docker -  command can be used to check.

You can test the container by doing the following:

First, make a test directory (just to keep things clean)
```bash
mkdir test
cd test
```

## Test your docker image

First, let's test to see if the java install in the docker image was successful

```bash
docker run -it -v $PWD:$PWD -w $PWD star STAR --version
```




