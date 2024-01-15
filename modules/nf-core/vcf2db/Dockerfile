FROM mambaorg/micromamba:1.4-focal

COPY ./vcf2db_env.yaml .

RUN micromamba update --file ./vcf2db_env.yaml --yes
RUN micromamba install --yes --name base --no-deps -c bioconda vcf2db==2020.02.24
RUN micromamba install --yes --name base --no-deps -c conda-forge coreutils
RUN micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"