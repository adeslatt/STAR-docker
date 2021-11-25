# Dockerfile for STAR
# https://github.com/alexdobin/STAR
FROM python:3

MAINTAINER Anne Deslattes Mays adeslat@scitechcon.org

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# install conda
ENV MINICONDA_VERSION py37_4.9.2
ENV CONDA_DIR /miniconda3

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash


# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# The conda bug with tbb - salmon: error while loading shared libraries: libtbb.so.2
# pandoc via conda was not working
RUN apt-get update && apt-get install -y libtbb2 pandoc-citeproc

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/star-2.7.2c/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name star-2.7.2c > star-2.7.2c.yml
