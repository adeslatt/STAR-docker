# Dockerfile for STAR
# https://github.com/alexdobin/STAR
FROM continuumio/miniconda3:4.8.2

RUN conda install -c bioconda star=2.7.2c
