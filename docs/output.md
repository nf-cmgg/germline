# nf-cmgg-germline: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

![Metro map of the pipeline workflow](images/nf-cmgg-germline_metro.png)

## VCFs

<!-- TODO update the output files for final release
-->

Following VCF files are currently used as output of the pipeline (subject to change):

1. Raw variant-called GVCFs and indexes in `<outdir>/individuals/<sample_name>`
2. The reblocked GVCFS in `<outdir>/individuals/<sample_name>`
3. The combined GVCFs and indexes in `<outdir>/families/<family_id>`
4. The genotyped or converted VCFs in `<outdir>/families/<family_id>`
5. The filtered VCFs (only for seqplorer mode) in `<outdir>/families/<family_id>`
6. The individual stat files in `<outdir>/families/<family_id>/reports`
7. The annotated VCFs (only for seqplorer mode) in `<outdir>/families/<family_id>`
8. The Gemini DB files (only for seqplorer mode) in `<outdir>/families/<family_id>`
9. The multiQC report in `<outdir>/multiqc_reports`

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
