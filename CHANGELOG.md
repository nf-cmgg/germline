# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev

- Added genomic data support
  - Added automatic BED file generation from the FASTA index
- Some minor improvements which could speed up pipeline execution (mainly scatter/gather improvements and additions)
- Added better support for the family and PED input (also support for samples that have neither input => will use the sample name as family name for a family containing this one sample)
- Refactored the code for better readability
- Added `dump` statements for better debugging
- Updated the minimum Nextflow version to 22.10.1
- Rewrote the post-processing subworkflow to a joint-genotyping workflow with GenomicsDBImport
- Added `dbsnp` input to Haplotypecaller
- Added `somalier` as a new subworkflow
- Improved `VCFanno`

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
