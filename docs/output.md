# nf-cmgg/germline: Output

## Introduction

This page describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level output directory (specified by `--outdir <DIR>`). This is an example output when the pipeline has been run for a WGS sample called `SAMPLE_1` and a WES sample called `SAMPLE_2` which form a family called `FAMILY_1`. The output consists of 4 directories: `yyyy-MM-dd_project_name`, `individuals`, `multiqc_reports` and `pipeline_info`. This run has only been run with `haplotypecaller`: (`--callers haplotypecaller`)

```bash
results/
├── Proband_12345
│   ├── NA24143_v2_0_0dev_2024_11_18
│   │   ├── NA24143.bed
│   │   ├── NA24143.haplotypecaller.bcftools_stats.txt
│   │   ├── NA24143.haplotypecaller.g.vcf.gz
│   │   ├── NA24143.haplotypecaller.g.vcf.gz.tbi
│   │   └── validation
│   │       └── haplotypecaller
│   │           ├── NA24143.fn.vcf.gz
│   │           ├── NA24143.fn.vcf.gz.tbi
│   │           ├── NA24143.fp.vcf.gz
│   │           ├── NA24143.fp.vcf.gz.tbi
│   │           ├── NA24143.non_snp.png
│   │           ├── NA24143.non_snp_roc.tsv.gz
│   │           ├── NA24143.non_snp.svg
│   │           ├── NA24143.phasing.txt
│   │           ├── NA24143.snp.png
│   │           ├── NA24143.snp_roc.tsv.gz
│   │           ├── NA24143.snp.svg
│   │           ├── NA24143.summary.txt
│   │           ├── NA24143.tp-baseline.vcf.gz
│   │           ├── NA24143.tp-baseline.vcf.gz.tbi
│   │           ├── NA24143.tp.vcf.gz
│   │           ├── NA24143.tp.vcf.gz.tbi
│   │           ├── NA24143.weighted.png
│   │           ├── NA24143.weighted_roc.tsv.gz
│   │           └── NA24143.weighted.svg
│   ├── NA24149_v2_0_0dev_2024_11_18
│   │   ├── NA24149.bed
│   │   ├── NA24149.haplotypecaller.bcftools_stats.txt
│   │   ├── NA24149.haplotypecaller.g.vcf.gz
│   │   ├── NA24149.haplotypecaller.g.vcf.gz.tbi
│   │   └── validation
│   │       └── haplotypecaller
│   │           ├── NA24149.fn.vcf.gz
│   │           ├── NA24149.fn.vcf.gz.tbi
│   │           ├── NA24149.fp.vcf.gz
│   │           ├── NA24149.fp.vcf.gz.tbi
│   │           ├── NA24149.non_snp.png
│   │           ├── NA24149.non_snp_roc.tsv.gz
│   │           ├── NA24149.non_snp.svg
│   │           ├── NA24149.phasing.txt
│   │           ├── NA24149.snp.png
│   │           ├── NA24149.snp_roc.tsv.gz
│   │           ├── NA24149.snp.svg
│   │           ├── NA24149.summary.txt
│   │           ├── NA24149.tp-baseline.vcf.gz
│   │           ├── NA24149.tp-baseline.vcf.gz.tbi
│   │           ├── NA24149.tp.vcf.gz
│   │           ├── NA24149.tp.vcf.gz.tbi
│   │           ├── NA24149.weighted.png
│   │           ├── NA24149.weighted_roc.tsv.gz
│   │           └── NA24149.weighted.svg
│   ├── NA24385_v2_0_0dev_2024_11_18
│   │   ├── NA24385.bed
│   │   ├── NA24385.haplotypecaller.bcftools_stats.txt
│   │   ├── NA24385.haplotypecaller.g.vcf.gz
│   │   ├── NA24385.haplotypecaller.g.vcf.gz.tbi
│   │   └── validation
│   │       └── haplotypecaller
│   │           ├── NA24385.fn.vcf.gz
│   │           ├── NA24385.fn.vcf.gz.tbi
│   │           ├── NA24385.fp.vcf.gz
│   │           ├── NA24385.fp.vcf.gz.tbi
│   │           ├── NA24385.non_snp.png
│   │           ├── NA24385.non_snp_roc.tsv.gz
│   │           ├── NA24385.non_snp.svg
│   │           ├── NA24385.phasing.txt
│   │           ├── NA24385.snp.png
│   │           ├── NA24385.snp_roc.tsv.gz
│   │           ├── NA24385.snp.svg
│   │           ├── NA24385.summary.txt
│   │           ├── NA24385.tp-baseline.vcf.gz
│   │           ├── NA24385.tp-baseline.vcf.gz.tbi
│   │           ├── NA24385.tp.vcf.gz
│   │           ├── NA24385.tp.vcf.gz.tbi
│   │           ├── NA24385.weighted.png
│   │           ├── NA24385.weighted_roc.tsv.gz
│   │           └── NA24385.weighted.svg
│   ├── output_v2_0_0dev_2024_11_18
│   │   ├── automap
│   │   │   └── haplotypecaller
│   │   │       ├── sample1
│   │   │       │   ├── sample1.HomRegions.cmgg_bio.tsv
│   │   │       │   ├── sample1.HomRegions.pdf
│   │   │       │   ├── sample1.HomRegions.strict.cmgg_bio.tsv
│   │   │       │   └── sample1.HomRegions.tsv
│   │   │       ├── sample2
│   │   │       │   ├── sample2.HomRegions.cmgg_bio.tsv
│   │   │       │   ├── sample2.HomRegions.pdf
│   │   │       │   ├── sample2.HomRegions.strict.cmgg_bio.tsv
│   │   │       │   └── sample2.HomRegions.tsv
│   │   │       └── sample3
│   │   │           ├── sample3.HomRegions.cmgg_bio.tsv
│   │   │           ├── sample3.HomRegions.pdf
│   │   │           ├── sample3.HomRegions.strict.cmgg_bio.tsv
│   │   │           └── sample3.HomRegions.tsv
│   │   ├── Proband_12345.haplotypecaller.bed
│   │   ├── Proband_12345.haplotypecaller.db
│   │   ├── Proband_12345.haplotypecaller.ped
│   │   ├── Proband_12345.haplotypecaller.vcf.gz
│   │   └── Proband_12345.haplotypecaller.vcf.gz.tbi
│   └── qc_v2_0_0dev_2024_11_18
│       ├── Proband_12345.haplotypecaller.bcftools_stats.txt
│       └── Proband_12345.haplotypecaller.html
└── v2_0_0dev_2024_11_18
    ├── execution_report_2024-11-18_14-03-56.html
    ├── execution_timeline_2024-11-18_14-03-56.html
    ├── execution_trace_2024-11-18_14-03-56.html
    ├── multiqc_report.html
    ├── params_2024-11-18_14-04-11.json
    ├── pipeline_dag_2024-11-18_14-03-56.html
    ├── pipeline_software_mqc_versions.yml
    └── samplesheet.csv
```


```bash
results/
├── YYYY_MM_DD_project_name #(1)!
│   └── FAMILY_1 #(2)!
│       ├── FAMILY_1.bed #(3)!
│       ├── FAMILY_1.haplotypecaller.ped #(4)!
│       ├── FAMILY_1.haplotypecaller.vcf.gz #(5)!
│       ├── FAMILY_1.haplotypecaller.vcf.gz.tbi #(6)!
│       └── reports
│           ├── FAMILY_1.haplotypecaller.bcftools_stats.txt #(7)!
│           └── FAMILY_1.haplotypecaller.somalier.html #(8)!
├── multiqc
│   ├── multiqc_data/
│   └── multiqc_report.html #(9)!
├── SAMPLE_1 #(10)!
│   ├── SAMPLE_1.bed #(11)!
│   ├── SAMPLE_1.haplotypecaller.g.vcf.gz #(12)!
│   ├── SAMPLE_1.haplotypecaller.g.vcf.gz.tbi
│   └── reports
│       ├── SAMPLE_1.haplotypecaller.bcftools_stats.txt
│       ├── SAMPLE_1.mosdepth.global.dist.txt #(13)!
│       └── SAMPLE_1.mosdepth.summary.txt #(14)!
├── SAMPLE_2
│   ├── SAMPLE_2.bed
│   ├── SAMPLE_2.haplotypecaller.g.vcf.gz
│   ├── SAMPLE_2.haplotypecaller.g.vcf.gz.tbi
│   └── reports
│       ├── SAMPLE_2.haplotypecaller.bcftools_stats.txt
│       ├── SAMPLE_2.mosdepth.global.dist.txt
│       └── SAMPLE_2.mosdepth.summary.txt
├── pipeline_info/ #(15)!
└── samplesheet.csv #(16)!
```

1. This is the name of the main pipeline output. It contains the current date and the mnemonic name of the pipeline run by default. The date can be excluded with the `--skip_date_project` parameter and the name can be customized with the `--project <STRING>` parameter.

2. This directory contains all files for family `FAMILY_1`.

3. This is the BED file used to parallelize the joint-genotyping. It contains all regions where real variants have been found in all GVCFs in the family. The value of `--merge_distance` (default: `100000` base pairs) is used to pad the region so the BED file contains multiple bigger regions instead of tons of small regions.

4. The PED file detailing the relation between the different members of the family. This file will be inferred when no PED file has been given to this family.

5. The resulting VCF for this family. All desired post-processing has been applied on this file.

6. The index of the resulting VCF.

7. The statistics created with `bcftools stats` for the resulting VCF.

8. The results of `somalier relate`.

9. The report created with MultiQC. This contains all statistics generated with `bcftools stats`, Ensembl VEP and other tools.

10. The folder for `SAMPLE_1` containing temporary files that could be useful for re-analysis later.

11. This is the BED file used to parallelize the variant calling. It contains all regions that are callable in the input files based on the desired regions (WGS = the whole genome; WES = the regions specified in the `roi` BED file).

12. The GVCF file created with `haplotypecaller`. This can be used in later runs of the pipeline to skip variant calling for this sample. A major use case for this is to add a new member to a family without having to call all variants of already called members.

13. The global distribution of the coverage calculated by `mosdepth`.

14. The summary created by `mosdepth`.

15. The directory containing information of the pipeline run.

16. The samplesheet used for the pipeline run.

## Pipeline overview

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
