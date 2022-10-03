# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. This is an example output when the pipeline has been run with the test data provided in the [samplesheet](../assets/samplesheet.csv). The output consists of 4 directories: `families`, `individuals`, `multiqc_reports` and `pipeline_info`.

- The folder `families` contains the combined VCFs of every individual in the same family along with the quality reports generated from these files
  - Seqr mode: The unfiltered VCFs after the merge of the individual VCFs, also contains the indices of these VCFs
  - Seqplorer mode: The filtered VCFs and the annotated VCFs. Also contains a Gemini DB file of the annotated VCFs
- The folder `individuals` contains the GVCF and index of every individual
- The folder `multiqc_reports` contains all the MultiQC report files
- The folder `pipeline_info` contains reports on the execution of the pipeline

### Seqr mode

```bash
results/
├── families
│   └── Proband_12345
│       ├── Proband_12345.vcf.gz
│       ├── Proband_12345.vcf.gz.tbi
│       └── reports
│           ├── Proband_12345.bcftools_stats.txt
│           ├── Proband_12345.FILTER.summary
│           ├── Proband_12345.TsTv.count
│           └── Proband_12345.TsTv.qual
├── individuals
│   ├── NA12878K12_NVQ_034
│   │   ├── NA12878K12_NVQ_034.g.vcf.gz
│   │   └── NA12878K12_NVQ_034.g.vcf.gz.tbi
│   └── NA24385D2_NVQ_034
│       ├── NA24385D2_NVQ_034.g.vcf.gz
│       └── NA24385D2_NVQ_034.g.vcf.gz.tbi
├── multiqc_reports
│   ├── multiqc_data
│   ├── multiqc_plots
│   └── multiqc_report.html
└── pipeline_info
    ├── execution_report_2022-10-03_11-56-25.html
    ├── execution_timeline_2022-10-03_11-56-25.html
    ├── execution_trace_2022-10-03_11-56-25.txt
    └── pipeline_dag_2022-10-03_11-56-25.html
```

### Seqplorer mode

```bash
results/
├── families
│   └── Proband_12345
│       ├── Proband_12345.ann.vcf.gz
│       ├── Proband_12345.db
│       ├── Proband_12345_filtered_snps_indels.vcf.gz
│       └── reports
│           ├── Proband_12345.bcftools_stats.txt
│           ├── Proband_12345.FILTER.summary
│           ├── Proband_12345.TsTv.count
│           └── Proband_12345.TsTv.qual
├── individuals
│   ├── NA12878K12_NVQ_034
│   │   ├── NA12878K12_NVQ_034.g.vcf.gz
│   │   └── NA12878K12_NVQ_034.g.vcf.gz.tbi
│   └── NA24385D2_NVQ_034
│       ├── NA24385D2_NVQ_034.g.vcf.gz
│       └── NA24385D2_NVQ_034.g.vcf.gz.tbi
├── multiqc_reports
│   ├── multiqc_data
│   ├── multiqc_plots
│   └── multiqc_report.html
└── pipeline_info
    ├── execution_report_2022-10-03_11-51-54.html
    ├── execution_timeline_2022-10-03_11-51-54.html
    ├── execution_trace_2022-10-03_11-51-54.txt
    └── pipeline_dag_2022-10-03_11-51-54.html
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

![Metro map of the pipeline workflow](images/nf-cmgg-germline_metro.png)

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
