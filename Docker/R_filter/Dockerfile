FROM r-base


RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

COPY filter.vh20200112.R /home/
COPY filter.vh20200112_singleEnd.R /home/
COPY NC_001416.1.fasta /home/

RUN apt-get update
RUN apt-get -y install  bedtools \
                        build-essential \
                        samtools \
                        libbz2-dev  \
                        libcurl4-openssl-dev \
                        liblzma-dev \
                        libssl-dev \
			libz-dev \
			zlib1g-dev \
			libxml2-dev
                                                
## install R-packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('data.table',dependencies=TRUE)"

RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "BiocManager::install('GenomicAlignments')"
