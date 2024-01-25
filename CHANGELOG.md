# CenterForMedicalGeneticsGhent/nf-cmgg-germline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.4.2 - Vibrant Veurne - [January 25 2024]

### Fixes

1. Set the default ensembl VEP version to 105.0 instead of using dynamic container fetching

## v1.4.1 - Lively Leuven - [January 15 2024]

### New Features

1. Added the `--output_suffix` parameter to add a custom suffix to the basename of the output files.
2. Implemented files for the alphamissense plugin of VEP.
3. Added the `--only_pass` parameter to only output variants that have the `PASS` flag in the FILTER column. (This is only applied when `--filter` is also given)
4. Added the `--keep_alt_contigs` parameter. This will tell the pipeline to not filter out the alternate contigs, which will now be done by default.
5. Add dbsnp Ids to VCFs coming from vardict. This will be done automatically if a dbsnp VCF is given to the pipeline through the `--dbsnp` parameter.

### Improvements

1. Updated the seqplorer profile so that the output filenames are correct for easy import
2. Changed the separator in `--vcfanno_resources` to `;`
   instead of `,` to allow commas in glob patterns.
3. Removed the reheader step from the vardict subworkflow and added a simple sed substitution to the vardictjava module
4. `vcf2db` now uses a python 2 environment to increase it's stability

## v1.4.0 - Kingly Kortrijk - [December 6 2023]

### New Features

1. Added the `--callers` parameter to specify the variant caller to use. Currently only `haplotypecaller` and `vardict` are supported.
2. Added the `vardict` variant caller.
3. Added the `--vardict_min_af` parameter to specify the minimum allele frequency for `vardict`. This option is also available in the samplesheet as `vardict_min_af` to set it dynamically per sample.
4. Added the `--output_genomicsdb` option to specify whether a GenomicsDB should be outputted or not. This will be `true` when using `only_merge`.
5. Added `--normalize` options for decomposing and normalizing of variants after calling and genotyping.
6. Added `WGS`, `WES`, `SeqCap`, `HyperCap` and `seqplorer` profiles that can be used to set the default parameters for these types of runs.

### Improvements

1. Refactored the pipeline to accomodate future additions of variant callers and genotypers
2. Removed a lot of unnecessary bloat
3. Improved GenomicsDBImport (can now be multithreaded and runs a lot faster). This will make very big runs more possible.
4. Changed `coverage_fast` to `mosdepth_slow`, reversing the effect of the parameter. By default mosdepth will now be run with `--fast-mode`. This can be disabled using the new `mosdepth_slow` parameter.
5. Automatically merge the regions that are within 150 bps of eachother for the variant calling. This way it's ensured that indel calling happens correctly.

### Fixes

1. Fixed an issue with the outputting of the validation PNG files, now all three types of PNGs are outputted.
2. Fixed a small issue where VCFs without a sample created by the callers could not be used by `bcftools concat`, these files will now be filtered from the input of the command.
3. Removed the `--maxentscan` parameter because this file is automatically present in the container

## v1.3.0 - Happy Hasselt - [July 10 2023]

### New Features

1. Added the `--only_call` parameter. Specifying this parameter tells the pipeline to only do variant calling and skip all post-processing. This will only output the GVCFs and files created to help variant calling.
2. The samplesheet is now also in the output folder.
3. Added an option `--only_merge` to tell the pipeline to create genomicsdbs and stop running there
4. Get regions from the GVCF instead of CRAM for joint genotyping. This removes the need to supply a CRAM file when a GVCF file has been used as input.

### Improvements

1. Updated `nf-validation` to v0.2.1.
2. Updated the samtools/merge tool to the nf-core version. This increases the efficiency and disk space usage of the tool.

### Fixes

1. Fixed an error where the truth VCFs caused a join error when the same sample was given multiple times
2. Updated some outdated error messages

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
3. Fixed a breaking bug where mosdepth didn't output the callable regions (this makes v1.2.0 deprecated, please use v1.2.1 instead)

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
