FROM mambaorg/micromamba:1.4-focal

# Install mamba dependencies
RUN micromamba install -y --name base -c conda-forge -c bioconda -c defaults \
    git \
    nim \
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

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

ENV somalier_version=0.2.15

# Install somalier
RUN curl -SL https://github.com/brentp/somalier/releases/download/v$somalier_version/somalier -o somalier && \
    chmod +x somalier && \
    cp somalier /usr/local/bin/

RUN which somalier && \
    somalier -h
