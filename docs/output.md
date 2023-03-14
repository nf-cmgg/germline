# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. This is an example output when the pipeline has been run with the test data provided in the [samplesheet](../assets/samplesheet.csv). The output consists of 4 directories: `ready`, `individuals`, `multiqc_reports` and `pipeline_info`.

- The folder `ready` contains the joint-genotyped VCFs of every family along with the output from `bcftools stats` and the somalier plots generated from these files. It also contains a PED file (automatically generated if it wasn't supplied in the samplesheet) and the extended version of the PED file as a `.samples.tsv` file. It also contains the BED regions used in the scatter/gather process for joint-genotyping.
- The folder `individuals` contains the GVCF and index of every individual along with the output of `bcftools stats` generated from these files. It also contains the BED regions used in the scatter/gather process for variant calling.
- The folder `multiqc_reports` contains all the MultiQC report files
- The folder `pipeline_info` contains reports on the execution of the pipeline

### Seqr mode

```bash
results/
├── ready/
│   └── FAMILY_1/
│       ├── FAMILY_1.vcf.gz
│       ├── FAMILY_1.vcf.gz.tbi
|       ├── FAMILY_1.ped
|       ├── FAMILY_1.samples.tsv
|       ├── FAMILY_1.windows.bed
│       └── reports/
│           └── FAMILY_1.bcftools_stats.txt
├── individuals/
│   ├── SAMPLE_1/
│   │   ├── SAMPLE_1.g.bed
│   │   ├── SAMPLE_1.g.vcf.gz
│   │   ├── SAMPLE_1.g.vcf.gz.tbi
│   │   ├── validation/
|   │   └── reports/
│   │       └── SAMPLE_1.bcftools_stats.txt
│   └── SAMPLE_2/
│       ├── SAMPLE_2.g.bed
│       ├── SAMPLE_2.g.vcf.gz
│       ├── SAMPLE_2.g.vcf.gz.tbi
│       ├── validation/
|       └── reports/
│           └── SAMPLE_2.bcftools_stats.txt
├── multiqc_reports/
│   ├── multiqc_data
│   ├── multiqc_plots
│   └── multiqc_report.html
└── pipeline_info/
    ├── execution_report_2022-10-03_11-56-25.html
    ├── execution_timeline_2022-10-03_11-56-25.html
    ├── execution_trace_2022-10-03_11-56-25.txt
    └── pipeline_dag_2022-10-03_11-56-25.html
```

## Pipeline overview

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
