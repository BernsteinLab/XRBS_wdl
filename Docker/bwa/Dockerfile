FROM ubuntu:22.04

RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

RUN apt-get update && \
    apt-get -y install  bedtools \
                        build-essential \
                        git \
                        libbz2-dev  \
                        libcurl4-openssl-dev \
                        liblzma-dev \
                        libssl-dev \
			libz-dev \
			zlib1g-dev
                                                
# Install bwa
RUN cd /tmp && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa; make && \
    cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin && \
    rm -fr /tmp/bwa

ARG GIT_TAG=__specify_in_build_arg__    
RUN echo $GIT_TAG > /VERSION
