FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

# Install tools
RUN apt-get update && apt-get install --yes openjdk-7-jre-headless \
    unzip

# set location
ENV FASTQC_LOC http://www.bioinformatics.babraham.ac.uk/projects/fastqc

# Download and install FastQC; replace low java memory requirement in fastq perl code
WORKDIR /usr/local
ADD ${FASTQC_LOC}/fastqc_v0.11.5.zip /tmp/
RUN unzip /tmp/fastqc_v0.11.5.zip && sed -i 's/Xmx250m/Xmx2048m/' FastQC/fastqc && chmod 755 FastQC/fastqc
ENV PATH /usr/local/FastQC/:$PATH

# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* 
