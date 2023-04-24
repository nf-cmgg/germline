FROM mambaorg/micromamba:1.4-focal

ENV SOMALIER_VERSION=0.2.15
ENV PIPELINE_VERSION=1.2.0

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
    bcftools=1.17 \
    ensembl-vep=105.0 \
    gatk4=4.4.0.0 \
    rtg-tools=3.12.1 \
    sed=4.7 \
    grep=3.4 \
    vcf2db=2020.02.24 \
    vcfanno=0.3.3 \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN nextflow pull CenterForMedicalGeneticsGhent/nf-cmgg-germline -r ${PIPELINE_VERSION} && \
    nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-germline -r ${PIPELINE_VERSION} --help

USER root

# Install somalier
RUN curl -SL https://github.com/brentp/somalier/releases/download/v${SOMALIER_VERSION}/somalier -o somalier && \
    chmod +x somalier && \
    cp somalier /usr/local/bin/ && \
    somalier -h
