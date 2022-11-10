# Dockerfile for STAR
# https://github.com/alexdobin/STAR
FROM ubuntu
LABEL description="Base docker image for STAR starting with ubuntu"
ARG ENV_NAME="star"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get -y install wget && \
    apt-get -y install unzip && \
    rm -rf /var/lib/apt/lists/*

#
# Lets download the latest release of STAR which on this date i 2.7.10b
#
RUN cd /usr/bin && \
    wget https://github.com/alexdobin/STAR/releases/download/2.7.10b/STAR_2.7.10b.zip && \
    unzip STAR_2.7.10b.zip && \
    ln -s STAR_2.7.10b/Linux_x86_64_static/STAR . && \
    ln -s STAR_2.7.10b/Linux_x86_64_static/STARlong .
