# CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters

A nextflow pipeline for calling and annotating variants

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with samples, and a header row. See [usage docs](./usage.md).</small></details>| `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genome` | Reference genome build <details><summary>Help</summary><small>Requires a Genome Reference Consortium reference ID (e.g. GRCh38)</small></details>| `string` | GRCh38 |  |  |
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified. The path to the reference genome fasta.</small></details>| `string` |  | True |  |
| `fai` | Path to FASTA genome index file. | `string` |  |  |  |
| `dict` | Path to the sequence dictionary generated from the FASTA reference | `string` |  |  |  |
| `strtablefile` | Path to the STR table file generated from the FASTA reference | `string` |  |  |  |
| `sdf` | Path to the SDF folder generated from the reference FASTA file | `string` |  |  |  |
| `genomes_base` | Directory base for CMGG reference store (used when --genomes_ignore false is specified) | `string` | /references/ |  |  |
| `cmgg_config_base` | The base directory for the local config files | `string` | /conf/ |  | True |
| `genomes_ignore` | Do not load the local references from the path specified with --genomes_base | `boolean` |  |  | True |
| `igenomes_base` | Directory / URL base for iGenomes references. | `string` | s3://ngi-igenomes/igenomes |  | True |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` | True |  | True |

## Pipeline specific parameters

Parameters that define how the pipeline works

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `scatter_count` | The amount of scattering that should happen per sample. <details><summary>Help</summary><small>Increase this number to increase the pipeline run speed, but at the tradeoff of using more IO and disk space. This can differ from the actual scatter count in some cases (especially with smaller files).<br>This has an effect on HaplotypeCaller, GenomicsDBImport and GenotypeGVCFs.</small></details>| `integer` | 40 | True |  |
| `merge_distance` | The merge distance for genotype BED files <details><summary>Help</summary><small>Increase this parameter if GenomicsDBImport is running slow. This defines the maximum distance between intervals that should be merged. The less intervals GenomicsDBImport actually gets, the faster it will run.</small></details>| `integer` | 10000 |  |  |
| `dragstr` | Create DragSTR models to be used with HaplotypeCaller <details><summary>Help</summary><small>This currently is only able to run single-core per sample. Due to this, the process is very slow with only very small improvements to the analysis.</small></details>| `boolean` |  |  |  |
| `validate` | Validate the found variants <details><summary>Help</summary><small>This only validates individual sample GVCFs that have truth VCF supplied to them via the samplesheet (in row `truth_vcf`, with an optional index in the `truth_tbi` row)</small></details>| `boolean` |  |  |  |
| `filter` | Filter the found variants | `boolean` |  |  |  |
| `annotate` | Annotate the found variants | `boolean` |  |  |  |
| `add_ped` | Add PED INFO header lines to the final VCFs | `boolean` |  |  |  |
| `gemini` | Create a Gemini databases from the final VCFs | `boolean` |  |  |  |
| `mosdepth_slow` | Don't run mosdepth in fast-mode <details><summary>Help</summary><small>This is advised if you need exact coverage BED files as output</small></details>| `boolean` |  |  |  |
| `project` | The name of the project. <details><summary>Help</summary><small>This will be used to specify the final output files folder in the output directory.</small></details>| `string` |  |  |  |
| `roi` | Path to the default ROI (regions of interest) BED file to be used for WES analysis <details><summary>Help</summary><small>This will be used for all samples that do not have a specific ROI file supplied to them through the samplesheet. Don't supply an ROI file to run the analysis as WGS.</small></details>| `string` |  |  |  |
| `dbsnp` | Path to the dbSNP VCF file | `string` |  |  |  |
| `dbsnp_tbi` | Path to the index of the dbSNP VCF file | `string` |  |  |  |
| `somalier_sites` | Path to the VCF file with sites for Somalier to use | `string` | https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz |  |  |
| `only_call` | Only call the variants without doing any post-processing | `boolean` |  |  |  |
| `only_merge` | Only run the pipeline until the creation of the genomicsdbs and output them | `boolean` |  |  |  |
| `output_genomicsdb` | Output the genomicsDB together with the joint-genotyped VCF | `boolean` |  |  |  |
| `callers` | A comma delimited string of the available callers. Current options are: 'haplotypecaller' and 'vardict' | `string` | haplotypecaller |  |  |
| `vardict_min_af` | The minimum allele frequency for VarDict when no `vardict_min_af` is supplied in the samplesheet | `number` | 0.1 |  |  |
| `normalize` | Normalize the VCF after joint genotyping (will run on the decomposed VCF when --decompose is also used) | `boolean` |  |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  |  |
| `version` | Display version and exit. | `boolean` |  |  |  |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  |  |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  |  |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  |  |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `validationShowHiddenParams` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `validationLenientMode` | Lenient mode for parameter validation | `boolean` | False |  | True |
| `validationFailUnrecognisedParams` | Fail on unrecognised parameters | `boolean` |  |  | True |
| `validationSchemaIgnoreParams` | Comma-separated list of parameters to ignore when validating against the schema | `string` | genomes,test_data |  | True |

## Annotation parameters

Parameters to configure Ensembl VEP and VCFanno

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `vep_chunk_size` | The amount of sites per split VCF as input to VEP | `integer` | 50000 |  |  |
| `species` | The species of the samples <details><summary>Help</summary><small>Must be lower case and have underscores as spaces</small></details>| `string` | homo_sapiens |  |  |
| `vep_merged` | Specify if the VEP cache is a merged cache | `boolean` | True |  |  |
| `vep_cache` | The path to the VEP cache | `string` | None |  |  |
| `vep_dbnsfp` | Use the dbNSFP plugin with Ensembl VEP <details><summary>Help</summary><small>The '--dbnsfp' and '--dbnsfp_tbi' parameters need to be specified when using this parameter.</small></details>| `boolean` |  |  |  |
| `vep_spliceai` | Use the SpliceAI plugin with Ensembl VEP <details><summary>Help</summary><small>The '--spliceai_indel', '--spliceai_indel_tbi', '--spliceai_snv' and '--spliceai_snv_tbi' parameters need to be specified when using this parameter.</small></details>| `boolean` |  |  |  |
| `vep_spliceregion` | Use the SpliceRegion plugin with Ensembl VEP | `boolean` |  |  |  |
| `vep_mastermind` | Use the Mastermind plugin with Ensembl VEP <details><summary>Help</summary><small>The '--mastermind' and '--mastermind_tbi' parameters need to be specified when using this parameter.</small></details>| `boolean` |  |  |  |
| `vep_maxentscan` | Use the MaxEntScan plugin with Ensembl VEP <details><summary>Help</summary><small>The '--maxentscan' parameter need to be specified when using this parameter.</small></details>| `boolean` |  |  |  |
| `vep_eog` | Use the custom EOG annotation with Ensembl VEP <details><summary>Help</summary><small>The '--eog' and '--eog_tbi' parameters need to be specified when using this parameter.</small></details>| `boolean` |  |  |  |
| `vep_version` | The version of the VEP tool to be used | `string` | 105.0 |  |  |
| `vep_cache_version` | The version of the VEP cache to be used | `string` | 105 |  |  |
| `dbnsfp` | Path to the dbSNFP file | `string` |  |  |  |
| `dbnsfp_tbi` | Path to the index of the dbSNFP file | `string` |  |  |  |
| `spliceai_indel` | Path to the VCF containing indels for spliceAI | `string` |  |  |  |
| `spliceai_indel_tbi` | Path to the index of the VCF containing indels for spliceAI | `string` |  |  |  |
| `spliceai_snv` | Path to the VCF containing SNVs for spliceAI | `string` |  |  |  |
| `spliceai_snv_tbi` | Path to the index of the VCF containing SNVs for spliceAI | `string` |  |  |  |
| `mastermind` | Path to the VCF for Mastermind | `string` |  |  |  |
| `mastermind_tbi` | Path to the index of the VCF for Mastermind | `string` |  |  |  |
| `eog` | Path to the VCF containing EOG annotations | `string` |  |  |  |
| `eog_tbi` | Path to the index of the VCF containing EOG annotations | `string` |  |  |  |
| `vcfanno` | Run annotations with vcfanno | `boolean` |  |  |  |
| `vcfanno_config` | The path to the VCFanno config TOML | `string` |  |  |  |
| `vcfanno_lua` | The path to a Lua script to be used in VCFanno | `string` |  |  |  |
| `vcfanno_resources` | A comma-seperated list of resource files for VCFanno, please also supply their indices using this parameter | `string` |  |  |  |
