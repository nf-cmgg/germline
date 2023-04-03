FROM mambaorg/micromamba:1.4-focal

# Install mamba dependencies
RUN micromamba install -y --name base -c conda-forge -c bioconda -c defaults \
    mosdepth=0.3.3 \
    multiqc=1.14 \
    nextflow=22.10.6 \
    samtools=1.16.1 \
    tar=1.34 \
    tabix=1.11 \
    bedtools=2.30.0 \
    bcftools=1.16 \
    ensembl-vep=105.0 \
    gatk4=4.4.0.0 \
    rtg-tools=3.12.1 \
    sed=4.7 \
    grep=3.4 \
    vcf2db=2020.02.24 \
    vcfanno=0.3.3 \
    && micromamba clean --all --yes

USER root

RUN apt update -y && apt install git nim wget libz-dev liblzma-dev libbz2-dev libcurl4-openssl-dev -y

RUN cd /usr/bin && \
    wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -vxjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    make

# Install somalier
RUN cd / &&    \
    git clone -b master --depth 5 https://github.com/brentp/somalier.git

RUN cd /somalier &&  \
    nim c -d:danger -d:nsb_static -d:release -d:openmp -d:blas=openblas -d:lapack=openblas -o:/usr/bin/somalier src/somalier && \
    cp scripts/ancestry-labels-1kg.tsv / && \
    rm -rf /somalier && somalier --help

ENV somalier_ancestry_labels /ancestry_labels-1kg.tsv