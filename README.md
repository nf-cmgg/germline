# ![nf-cmgg-germline](docs/images/nf-cmgg-germline_logo_light.png#gh-light-mode-only) ![nf-cmgg-germline](docs/images/nf-cmgg-germline_logo_dark.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

## Introduction

**nf-cmgg-germline** is a bioinformatics best-practice analysis pipeline for A nextflow pipeline for calling and annotating variants. It uses HaplotypeCaller to call variants and EnsemblVEP to annotate the called variants. By supplying the `--output_mode <seqr|seqplorer>` you can choose for which platform the VCFs should be created.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Check the input CSV file (This checks the CSV format, headers and files passed throught the CSV)
2. (Only if `--fasta_fai FILE` is not used) Create a FASTA index file from the FASTA reference ([Samtools Faidx](http://www.htslib.org/doc/samtools-faidx.html))
3. (Only if `--dict FILE` is not used) Create a sequence dictionary file from the FASTA reference ([GATK CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-))
4. (Only if `--strtablefile FILE` is not used and `--use_dragstr_model true` is given) Create an STR table file from the FASTA reference ([GATK ComposeSTRTableFile](https://gatk.broadinstitute.org/hc/en-us/articles/4405451249819-ComposeSTRTableFile))
5. (Only if `--scatter_count INT` is bigger than 1) Split the BED files for parallellization of HaplotypeCaller ([Bedtools Split](https://bedtools.readthedocs.io/en/latest/content/overview.html))
6. (Only if `--use_dragstr_model true` is given) Create a DRAGstr model for each sample provided ([GATK CalibrateDragstrModel](https://gatk.broadinstitute.org/hc/en-us/articles/360057441571-CalibrateDragstrModel-BETA-))
7. Call the variants for each sample. The BED files are used to parallellize this process if `--scatter_count INT` is bigger than 1 ([GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller))
8. (Only if `--scatter_count INT` is bigger than 1) Concatenate the GVCF files for each sample ([Bcftools Concat](https://samtools.github.io/bcftools/bcftools.html#concat))
9. Reblock the GVCFs created with HaplotypeCaller ([GATK ReblockGVCF](https://gatk.broadinstitute.org/hc/en-us/articles/4405443600667-ReblockGVCF))
10. Combine the reblocked GVCFs of the same family ([GATK CombineGVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs))
11. Genotype the combined GVCFs ([GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs))
12. Create a VCF header with pedigree information extracted from the PED file ([RTG pedfilter](https://www.animalgenome.org/bioinfo/resources/manuals/RTGOperationsManual.pdf) => chapter 2.5.18)
13. Merge the pedigree header with the header of the genotyped VCF (This is done using a custom-written python [script](bin/merge_vcf_headers.py))
14. (Only if `--output_mode <seqr|seqplorer>` is set to `seqplorer`) Filter the VCF to be the correct format for Seqplorer ([bcftools filter](http://samtools.github.io/bcftools/bcftools.html#filter))
15. Perform a quality control on the VCFs ([bcftools stats](http://samtools.github.io/bcftools/bcftools.html#stats) and [vcftools](http://vcftools.sourceforge.net/man_latest.html))
16. (Only if `--output_mode <seqr|seqplorer>` is set to `seqplorer`) Annotate the genotyped VCF files ([Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html))
17. (Only if `--output_mode <seqr|seqplorer>` is set to `seqplorer`) Transform the VCF to a Gemini-compatible database file for Seqplorer compatibility ([vcf2db](https://github.com/quinlan-lab/vcf2db))
18. Run [MultiQC](https://multiqc.info/) on all quality control files

![metro graph](docs/images/nf-cmgg-germline_metro.png)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-germline -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```console
   nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-germline --input <INPUT_CSV> --outdir <OUTDIR> --genome GRCh38 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --fasta <PATH/TO/FASTA>
   ```

An overview of the parameters for this pipeline can be viewed using:

```
nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-germline --help
```

## Credits

nf-cmgg-germline was originally written by @nvnieuwk.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
