# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level output directory (specified by `--outdir <DIR>`). This is an example output when the pipeline has been run for a WGS sample called `SAMPLE_1` and a WES sample called `SAMPLE_2` which form a family called `FAMILY_1`. The output consists of 4 directories: `yyyy-MM-dd_project_name`, `individuals`, `multiqc_reports` and `pipeline_info`.

- The folder `yyyy-MM-dd_project_name` contains the final outputs (either joint-genotyped VCFs or single-sample VCFs dependent on what callers have been used). Alongside the VCFs, a generated PED file, BED file for joint-genotyped samples and gemini database (`--gemini`) are also present. The `yyyy-MM-dd` part of the folder will be the date of the analysis and `project_name` can be specified by `--project <STRING>`. If this parameter isn't specified, the mnemonic name of the nextflow run will be used instead.
- The sample folders contain the GVCF and index of every individual along with the output of `bcftools stats` generated from these files and the mosdepth output. It also contains the BED regions used in the scatter/gather process for variant calling. If a validation has been performed all results will also be in the subfolder `validation`.
- The folder `multiqc_reports` contains all the MultiQC report files
- The folder `pipeline_info` contains reports on the execution of the pipeline

```bash
results/
├── yyyy-MM-dd_project_name/
│   ├── FAMILY_1/
│   │   ├── FAMILY_1.caller.vcf.gz
│   │   ├── FAMILY_1.caller.vcf.gz.tbi
│   │   ├── FAMILY_1.caller.ped
│   │   ├── FAMILY_1.caller.db
│   │   ├── FAMILY_1.bed
│   │   └── reports/
│   │       ├── FAMILY_1.caller.bcftools_stats.txt
│   │       └── FAMILY_1.caller.somalier.html
│   ├── SAMPLE_1/
│   │   ├── SAMPLE_1.caller.vcf.gz
│   │   ├── SAMPLE_1.caller.vcf.gz.tbi
│   │   ├── SAMPLE_1.caller.ped
│   │   ├── SAMPLE_1.caller.db
│   │   └── reports/
│   │       ├── SAMPLE_1.caller.bcftools_stats.txt
│   │       └── SAMPLE_1.caller.somalier.html
│   └── SAMPLE_2/
│       ├── SAMPLE_2.caller.vcf.gz
│       ├── SAMPLE_2.caller.vcf.gz.tbi
│       ├── SAMPLE_2.caller.ped
│       ├── SAMPLE_2.caller.db
│       └── reports/
│           ├── SAMPLE_2.caller.bcftools_stats.txt
│           └── SAMPLE_2.caller.somalier.html
├── SAMPLE_1/
│   ├── SAMPLE_1.quantized.bed.gz
│   ├── SAMPLE_1.quantized.bed.gz.csi
│   ├── SAMPLE_1.mosdepth.global.dist.txt
│   ├── SAMPLE_1.summary.txt
│   ├── SAMPLE_1.caller.g.vcf.gz
│   ├── SAMPLE_1.caller.g.vcf.gz.tbi
│   ├── validation/
│   └── reports/
│       └── SAMPLE_1.bcftools_stats.txt
├── SAMPLE_2/
│   ├── SAMPLE_2.quantized.bed.gz
│   ├── SAMPLE_2.quantized.bed.gz.csi
│   ├── SAMPLE_2.regions.bed.gz
│   ├── SAMPLE_2.regions.bed.gz.csi
│   ├── SAMPLE_2.mosdepth.global.dist.txt
│   ├── SAMPLE_2.mosdepth.region.dist.txt
│   ├── SAMPLE_2.summary.txt
│   ├── SAMPLE_2.bed
│   ├── SAMPLE_2.caller.g.vcf.gz
│   ├── SAMPLE_2.caller.g.vcf.gz.tbi
│   ├── validation/
│   └── reports/
│       └── SAMPLE_2.bcftools_stats.txt
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
