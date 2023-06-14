# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.3dev

### New Features

Added the `--only_call` parameter. Specifying this parameter tells the pipeline to only do variant calling and skip all post-processing. This will only output the GVCFs and files created to help variant calling.

### Improvements

Updated `nf-validation` to v0.2.1.

## v1.2.2 - Benign Brussels - [June 12 2023]

### Improvements

Changed the output directory structure to be more bcbio like

## v1.2.1 - Balanced Brussels - [June 5 2023]

### New Features

1. Added support for the [`nf-validation`](https://github.com/nextflow-io/nf-validation/tree/master) plugin.
2. Haplotypecaller dragen mode will be automatically disabled when not using a dragstr model.

### Bug fixes

1. Removed bedtools/jaccard
2. Fixed some patterns in the parameter JSON schema (since they are actually used now)
3. Fixed a breaking bug where mosdepth didn't output the callable regions (this makes v1.2.0 deprecated, please use v1.3.0 instead)

### Improvements

1. Genomicsdbs aren't scattered now, this increases the precision of the analyis by almost 3% at the cost of a bit longer runtimes
2. Actually do the validation on the output VCFs now instead the freshly called GVCFs
3. Improved the efficiency of the VEP run by scattering more efficiently on the amount of variants instead of the chromosomes

## v1.2.0 - Brave Brussels - [May 5 2023]

### New Features

1. Added a `--coverage_fast <true/false>` flag which can be used to run mosdepth in fast mode. This flag will also make sure that only the quantized bed from mosdepth is present in the output directory for each WGS individual, otherwise it will output everything
2. Added the possibility to give GVCF files as inputs and immediately go to the joint-genotyping. This is especially useful for the cases where several samples should be combined. This way the variant calling doesn't need to be re-run. Beware though that a CRAM file should still be given to generate the BED files used for the scatter/gathering. The new header names are `gvcf` and `tbi` where `gvcf` is used to give the GVCF and `tbi` is used to give its index.
3. Added `bedtools jaccard` to the validation.
4. Added a Dockerfile which creates an image that is able to run a full pipeline run inside of it.
5. Added better documentation

### Improvements

1. Updated the scattering again: it now follows this workflow:
   - Sort and merge overlapping intervals of given ROI BED files (WES only)
   - Create a BED file with callable regions using mosdepth
   - Intersect the callable regions BED with the ROI BED (WES only)
   - Split the resulting BED file (or the callable regions BED for WGS) into evenly sized BED files (amount is specified with `--scatter_count`)
   - Run HaplotypeCaller in parallel using these regions
   - Merge and sort the BED files of all individuals in a family
   - Split the merged BED file into evenly sized BED files (amount is specified with `--scatter_count` times the family size)
   - Run GenomicsDBImport and GenotypeGVCFs in parallel using these regions
2. Updated the resource requirements of GenomicsDBImport and GenotypeGVCFs to be more efficient (and more cluster friendly)
3. Removed ReblockGVCFs (this wasn't worth it and we save the raw GVCFs)
4. Added `--merge_distance <integer>` to decrease the amount of intervals passed to genomicsdbimport. Increase this value if GenomicsDBImport is running slow.
5. Renamed `--use_dragstr_model` to `--dragstr`.

### Bug fixes

1. Fixed a warning showing up when running with `--dragstr false`
2. Add `--infer` flag to `somalier relate` when no PED file is given

## v1.1.2 - Groovy Ghent - [Mar 21 2023]

### New features

1. Added a parameter for setting the splitting depth threshold `--split_threshold FLOAT`

### Changes

1. Change the default splitting threshold to 0.2 instead of 0.3

## v1.1.1 - Golden Ghent - [Mar 20 2023]

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
