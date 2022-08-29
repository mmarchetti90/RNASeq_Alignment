FROM ubuntu:latest

ARG STAR_VERSION=2.7.10a
ARG SAMTOOLS_VERSION=1.16

## Basic resources
RUN apt-get -q update && \
    apt-get -q -y --no-install-recommends install \
        bzip2 \
        default-jre \
        gcc \
        g++ \
	libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        make \
        wget \
        zlib1g-dev && \
    apt-get clean

## FastQC
RUN apt-get -q -y install fastqc && \
    apt-get clean

## STAR
RUN cd /tmp && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && \
    tar -xzf ${STAR_VERSION}.tar.gz -C /home && \
    cd /home/STAR-${STAR_VERSION}/source && \
    make STAR && \
    rm /tmp/${STAR_VERSION}.tar.gz

## samtools
RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /home && \
    cd /home/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install && \
    rm /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2

## Adding STAR and samtools to PATH
ENV PATH $PATH:/home/STAR-${STAR_VERSION}/source:/home/samtools-${SAMTOOLS_VERSION}

## Set workdir to /home/
WORKDIR /home/

## Launch R automatically
CMD ["/bin/bash"]