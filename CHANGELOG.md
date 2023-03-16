# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.1dev

### Changes

1. Set the default of `--validate` to `false`

### Bug fixes

1. Fixed a bug with ensembl VEP. Filenames of the alt contigs should now have a `_alt` suffix instead of all alt contigs.
2. Added file-exist check to the `sdf` file
3. Fixed the scattering when using alt contigs

## v1.1.0 - Glorious Ghent - [Mar 14 2023]

### New Features

1. BED file input is new optional (The regions are created from the FASTA index). Providing a BED file is still preffered for the most optimal runs.
2. Added support for samples that aren't part of a family. Just leave the `ped` and `family_id` input fields in the samplesheet empty for a sample to be treated like this. This sample will go through exactly the same workflow but will be emitted as a single-sample VCF.
3. Added `dump` functionality to lots of channels.
4. Added the `dbsnp` option to `GATK HaplotypeCaller`. use `--dbsnp` and `--dbsnp_tbi` to supply these VCFs.
5. Added the `vcf_extract_somalier` subworkflow to the pipeline. This also creates PED files inferred from the input multi-sample VCF.
6. Added a validation subworkflow. All files that have a VCF in the `truth_vcf` column of the input samplesheet will be validated against this VCF. This can be turned off by supplying the `--validate false` flag to the pipeline run.

### Improvements

1. Improved the scatter/gather logic. This is now done with `goleft indexsplit` to define chunks of even coverage. The genotyping scattering now happens with `bedtools makewindows`. This creates chunks of even regions from the merged BED files for the family. By passing a padding of about 20 bps to the genotype tools, we make sure all variants on the edges of these regions are also genotyped. Duplicates are removed later when running `bcftools concat`
2. Refactored a lot of the code to maintain the same style over the whole pipeline.
3. Updated the minimum Nextlow version to `22.10.5` to make sure S3 staging works perfectly.
4. The `post_processing` subworklow has been renamed to the better suiting `joint_genotyping` subworkflow. `reblockgvcf` has been moved to `germline_variant_calling` and the `filter` and `reheadering` has been moved to the main workflow.
5. Merging VCFs of the same family now happens with `GATK GenomicsDBImport` instead of `GATK MergeGVCFs` or `bcftools merge`. This gives more reliable results.
6. Improved the handling of `vcfanno`
7. The PED headers can now be added to all the output VCFs that are part of a family instead of only those that were given a PED file as input. The PED file used is created using `somalier relate`. This feature can be turned on using the `--add_ped true` argument. This doesn't happen by default.

### Bug fixes

1. Fixed some issues when both the `ped` and `family_id` were given for a sample.
2. Fixed the PED input for `rtgtools_pedfilter` (`-9` isn't recognized as unknown by the tool. Now these will be automatically converted to `0` before this tool)
3. Fixed issues with the DBsnp index not being created correctly
4. Fixed wrongly formed joins and added checks for mismatches and duplicates

## v1.0.1 - Happy Hollebeke - [Oct 7 2022]

### Changes

- Upgraded to `nf-core` v2.6 template

### Fixes

- Fixed the `ensemblvep` version (was 104.1 before and is now 105.0)
- Updated the label of `gatk4/calibratedragstrmodel` to `process_high` to match the requirements for bigger inputs

## v1.0.0 - Beautiful Bruges - [Oct 3 2022]

### Added

- Full release of the pipeline

## v1.0dev - [May 31 2022]

Initial release of CenterForMedicalGeneticsGhent/nf-cmgg-germline, created with the [nf-core](https://nf-co.re/) template.
