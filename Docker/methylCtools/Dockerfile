FROM python:2.7

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
                        zlib1g-dev \
                        samtools \
                        tabix

RUN pip install pysam==0.13

# Install bwa
RUN git clone https://github.com/hovestadt/methylCtools.git 

RUN chmod +x methylCtools/methylCtools

ARG GIT_TAG=__specify_in_build_arg__    
RUN echo $GIT_TAG > /VERSION
