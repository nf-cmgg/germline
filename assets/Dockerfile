FROM mambaorg/micromamba:1.4-focal

ENV SOMALIER_VERSION=0.2.15

# Install mamba dependencies
RUN micromamba install -y --name base -c conda-forge -c bioconda -c defaults \
    git \
    nim \
    gcc \
    snappy \
    nextflow=22.10.6 \
    python==3.11 \
    mosdepth=0.3.3 \
    multiqc=1.14 \
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
    vcfanno=0.3.3 \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

# Install somalier
RUN curl -SL https://github.com/brentp/somalier/releases/download/v${SOMALIER_VERSION}/somalier -o somalier && \
    chmod +x somalier && \
    cp somalier /usr/local/bin/ && \
    somalier -h

# Install vcf2db
RUN git clone https://github.com/quinlan-lab/vcf2db && \
    cd vcf2db && \
    pip install -r requirements.txt && \
    pip install SQLAlchemy==1.4.4 && \
    chmod a+x vcf2db.py && \
    cp vcf2db.py /usr/local/bin/ && \
    vcf2db.py -h
