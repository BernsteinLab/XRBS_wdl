FROM ubuntu:22.04

RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

RUN apt-get update && \
    apt-get -y install  build-essential \
                        libbz2-dev  \
                        libcurl4-openssl-dev \
                        liblzma-dev \
                        libssl-dev \
                        libz-dev \
                        zlib1g-dev \
                        cutadapt \
                        curl

RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz

ARG GIT_TAG=__specify_in_build_arg__    
RUN echo $GIT_TAG > /VERSION
