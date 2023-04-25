# CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters

- [1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > input_output_options`](#allOf_i0)
  - [1.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Input/output options > input`](#allOf_i0_input)
  - [1.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Input/output options > outdir`](#allOf_i0_outdir)
  - [1.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Input/output options > email`](#allOf_i0_email)
- [2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > reference_genome_options`](#allOf_i1)
  - [2.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > genome`](#allOf_i1_genome)
  - [2.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > fasta`](#allOf_i1_fasta)
  - [2.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > fai`](#allOf_i1_fai)
  - [2.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > dict`](#allOf_i1_dict)
  - [2.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > strtablefile`](#allOf_i1_strtablefile)
  - [2.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > sdf`](#allOf_i1_sdf)
  - [2.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > genomes_base`](#allOf_i1_genomes_base)
  - [2.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > cmgg_config_base`](#allOf_i1_cmgg_config_base)
  - [2.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > genomes_ignore`](#allOf_i1_genomes_ignore)
  - [2.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > igenomes_base`](#allOf_i1_igenomes_base)
  - [2.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > igenomes_ignore`](#allOf_i1_igenomes_ignore)
- [3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > pipeline_specific_parameters`](#allOf_i2)
  - [3.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > scatter_count`](#allOf_i2_scatter_count)
  - [3.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > merge_distance`](#allOf_i2_merge_distance)
  - [3.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > dragstr`](#allOf_i2_dragstr)
  - [3.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > validate`](#allOf_i2_validate)
  - [3.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > filter`](#allOf_i2_filter)
  - [3.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > annotate`](#allOf_i2_annotate)
  - [3.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > add_ped`](#allOf_i2_add_ped)
  - [3.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > gemini`](#allOf_i2_gemini)
  - [3.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > coverage_fast`](#allOf_i2_coverage_fast)
  - [3.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > roi`](#allOf_i2_roi)
  - [3.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > dbsnp`](#allOf_i2_dbsnp)
  - [3.12. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > dbsnp_tbi`](#allOf_i2_dbsnp_tbi)
  - [3.13. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > somalier_sites`](#allOf_i2_somalier_sites)
- [4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > institutional_config_options`](#allOf_i3)
  - [4.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > custom_config_version`](#allOf_i3_custom_config_version)
  - [4.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > custom_config_base`](#allOf_i3_custom_config_base)
  - [4.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_name`](#allOf_i3_config_profile_name)
  - [4.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_description`](#allOf_i3_config_profile_description)
  - [4.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_contact`](#allOf_i3_config_profile_contact)
  - [4.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_url`](#allOf_i3_config_profile_url)
- [5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > max_job_request_options`](#allOf_i4)
  - [5.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Max job request options > max_cpus`](#allOf_i4_max_cpus)
  - [5.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Max job request options > max_memory`](#allOf_i4_max_memory)
  - [5.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Max job request options > max_time`](#allOf_i4_max_time)
- [6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > generic_options`](#allOf_i5)
  - [6.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > help`](#allOf_i5_help)
  - [6.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > version`](#allOf_i5_version)
  - [6.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > publish_dir_mode`](#allOf_i5_publish_dir_mode)
  - [6.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > email_on_fail`](#allOf_i5_email_on_fail)
  - [6.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > plaintext_email`](#allOf_i5_plaintext_email)
  - [6.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > max_multiqc_email_size`](#allOf_i5_max_multiqc_email_size)
  - [6.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > monochrome_logs`](#allOf_i5_monochrome_logs)
  - [6.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > hook_url`](#allOf_i5_hook_url)
  - [6.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_title`](#allOf_i5_multiqc_title)
  - [6.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_config`](#allOf_i5_multiqc_config)
  - [6.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_logo`](#allOf_i5_multiqc_logo)
  - [6.12. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_methods_description`](#allOf_i5_multiqc_methods_description)
  - [6.13. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > tracedir`](#allOf_i5_tracedir)
  - [6.14. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > validate_params`](#allOf_i5_validate_params)
  - [6.15. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > show_hidden_params`](#allOf_i5_show_hidden_params)
- [7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > annotation_parameters`](#allOf_i6)
  - [7.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > species`](#allOf_i6_species)
  - [7.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_merged`](#allOf_i6_vep_merged)
  - [7.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_cache`](#allOf_i6_vep_cache)
  - [7.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_dbnsfp`](#allOf_i6_vep_dbnsfp)
  - [7.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_spliceai`](#allOf_i6_vep_spliceai)
  - [7.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_spliceregion`](#allOf_i6_vep_spliceregion)
  - [7.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_mastermind`](#allOf_i6_vep_mastermind)
  - [7.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_maxentscan`](#allOf_i6_vep_maxentscan)
  - [7.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_eog`](#allOf_i6_vep_eog)
  - [7.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_version`](#allOf_i6_vep_version)
  - [7.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_cache_version`](#allOf_i6_vep_cache_version)
  - [7.12. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > dbnsfp`](#allOf_i6_dbnsfp)
  - [7.13. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > dbnsfp_tbi`](#allOf_i6_dbnsfp_tbi)
  - [7.14. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_indel`](#allOf_i6_spliceai_indel)
  - [7.15. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_indel_tbi`](#allOf_i6_spliceai_indel_tbi)
  - [7.16. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_snv`](#allOf_i6_spliceai_snv)
  - [7.17. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_snv_tbi`](#allOf_i6_spliceai_snv_tbi)
  - [7.18. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > mastermind`](#allOf_i6_mastermind)
  - [7.19. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > mastermind_tbi`](#allOf_i6_mastermind_tbi)
  - [7.20. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > maxentscan`](#allOf_i6_maxentscan)
  - [7.21. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > eog`](#allOf_i6_eog)
  - [7.22. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > eog_tbi`](#allOf_i6_eog_tbi)
  - [7.23. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno`](#allOf_i6_vcfanno)
  - [7.24. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno_config`](#allOf_i6_vcfanno_config)
  - [7.25. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno_lua`](#allOf_i6_vcfanno_lua)
  - [7.26. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno_resources`](#allOf_i6_vcfanno_resources)

**Title:** CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `combining`                                                               |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |

**Description:** A nextflow pipeline for calling and annotating variants

| All of(Requirement)                       |
| ----------------------------------------- |
| [input_output_options](#allOf_i0)         |
| [reference_genome_options](#allOf_i1)     |
| [pipeline_specific_parameters](#allOf_i2) |
| [institutional_config_options](#allOf_i3) |
| [max_job_request_options](#allOf_i4)      |
| [generic_options](#allOf_i5)              |
| [annotation_parameters](#allOf_i6)        |

## <a name="allOf_i0"></a>1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > input_output_options`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Defined in**            | #/definitions/input_output_options                                        |

**Description:** Define where the pipeline should find input data and save output data.

| Property                      | Pattern | Type   | Deprecated | Definition | Title/Description                                                                                                        |
| ----------------------------- | ------- | ------ | ---------- | ---------- | ------------------------------------------------------------------------------------------------------------------------ |
| + [input](#allOf_i0_input )   | No      | string | No         | -          | Path to comma-separated file containing information about the samples in the experiment.                                 |
| + [outdir](#allOf_i0_outdir ) | No      | string | No         | -          | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. |
| - [email](#allOf_i0_email )   | No      | string | No         | -          | Email address for completion summary.                                                                                    |

### <a name="allOf_i0_input"></a>1.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Input/output options > input`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | Yes         |
| **Format**   | `file-path` |

**Description:** Path to comma-separated file containing information about the samples in the experiment.

| Restrictions                      |                                                                                                                      |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.(csv\|tsv\|yaml\|yml)$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.%28csv%7Ctsv%7Cyaml%7Cyml%29%24) |

### <a name="allOf_i0_outdir"></a>1.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Input/output options > outdir`

|              |                  |
| ------------ | ---------------- |
| **Type**     | `string`         |
| **Required** | Yes              |
| **Format**   | `directory-path` |

**Description:** The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

### <a name="allOf_i0_email"></a>1.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Input/output options > email`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Email address for completion summary.

| Restrictions                      |                                                                                                                                                                                                                   |
| --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$``` [Test](https://regex101.com/?regex=%5E%28%5Ba-zA-Z0-9_%5C-%5C.%5D%2B%29%40%28%5Ba-zA-Z0-9_%5C-%5C.%5D%2B%29%5C.%28%5Ba-zA-Z%5D%7B2%2C5%7D%29%24) |

## <a name="allOf_i1"></a>2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > reference_genome_options`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Defined in**            | #/definitions/reference_genome_options                                    |

**Description:** Reference genome related files and options required for the workflow.

| Property                                          | Pattern | Type    | Deprecated | Definition | Title/Description                                                                       |
| ------------------------------------------------- | ------- | ------- | ---------- | ---------- | --------------------------------------------------------------------------------------- |
| - [genome](#allOf_i1_genome )                     | No      | string  | No         | -          | Reference genome build                                                                  |
| + [fasta](#allOf_i1_fasta )                       | No      | string  | No         | -          | Path to FASTA genome file.                                                              |
| - [fai](#allOf_i1_fai )                           | No      | string  | No         | -          | Path to FASTA genome index file.                                                        |
| - [dict](#allOf_i1_dict )                         | No      | string  | No         | -          | Path to the sequence dictionary generated from the FASTA reference                      |
| - [strtablefile](#allOf_i1_strtablefile )         | No      | string  | No         | -          | Path to the STR table file generated from the FASTA reference                           |
| - [sdf](#allOf_i1_sdf )                           | No      | string  | No         | -          | Path to the SDF folder generated from the reference FASTA file                          |
| - [genomes_base](#allOf_i1_genomes_base )         | No      | string  | No         | -          | Directory base for CMGG reference store (used when --genomes_ignore false is specified) |
| - [cmgg_config_base](#allOf_i1_cmgg_config_base ) | No      | string  | No         | -          | The base directory for the local config files                                           |
| - [genomes_ignore](#allOf_i1_genomes_ignore )     | No      | boolean | No         | -          | Do not load the local references from the path specified with --genomes_base            |
| - [igenomes_base](#allOf_i1_igenomes_base )       | No      | string  | No         | -          | Directory / URL base for iGenomes references.                                           |
| - [igenomes_ignore](#allOf_i1_igenomes_ignore )   | No      | boolean | No         | -          | Do not load the iGenomes reference config.                                              |

### <a name="allOf_i1_genome"></a>2.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > genome`

|              |            |
| ------------ | ---------- |
| **Type**     | `string`   |
| **Required** | No         |
| **Default**  | `"GRCh38"` |

**Description:** Reference genome build

### <a name="allOf_i1_fasta"></a>2.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > fasta`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | Yes         |
| **Format**   | `file-path` |

**Description:** Path to FASTA genome file.

| Restrictions                      |                                                                                                 |
| --------------------------------- | ----------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.fn?a(sta)?$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.fn%3Fa%28sta%29%3F%24) |

### <a name="allOf_i1_fai"></a>2.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > fai`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to FASTA genome index file.

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.fai$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.fai%24) |

### <a name="allOf_i1_dict"></a>2.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > dict`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the sequence dictionary generated from the FASTA reference

| Restrictions                      |                                                                             |
| --------------------------------- | --------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.dict$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.dict%24) |

### <a name="allOf_i1_strtablefile"></a>2.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > strtablefile`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |
| **Format**   | `path`   |

**Description:** Path to the STR table file generated from the FASTA reference

### <a name="allOf_i1_sdf"></a>2.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > sdf`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |
| **Format**   | `path`   |

**Description:** Path to the SDF folder generated from the reference FASTA file

### <a name="allOf_i1_genomes_base"></a>2.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > genomes_base`

|              |                  |
| ------------ | ---------------- |
| **Type**     | `string`         |
| **Required** | No               |
| **Format**   | `directory-path` |
| **Default**  | `"/references/"` |

**Description:** Directory base for CMGG reference store (used when --genomes_ignore false is specified)

### <a name="allOf_i1_cmgg_config_base"></a>2.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > cmgg_config_base`

|              |            |
| ------------ | ---------- |
| **Type**     | `string`   |
| **Required** | No         |
| **Default**  | `"/conf/"` |

**Description:** The base directory for the local config files

### <a name="allOf_i1_genomes_ignore"></a>2.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > genomes_ignore`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Do not load the local references from the path specified with --genomes_base

### <a name="allOf_i1_igenomes_base"></a>2.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > igenomes_base`

|              |                                |
| ------------ | ------------------------------ |
| **Type**     | `string`                       |
| **Required** | No                             |
| **Format**   | `directory-path`               |
| **Default**  | `"s3://ngi-igenomes/igenomes"` |

**Description:** Directory / URL base for iGenomes references.

### <a name="allOf_i1_igenomes_ignore"></a>2.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Reference genome options > igenomes_ignore`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |
| **Default**  | `true`    |

**Description:** Do not load the iGenomes reference config.

## <a name="allOf_i2"></a>3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > pipeline_specific_parameters`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Default**               | `""`                                                                      |
| **Defined in**            | #/definitions/pipeline_specific_parameters                                |

**Description:** Parameters that define how the pipeline works

| Property                                      | Pattern | Type    | Deprecated | Definition | Title/Description                                                                  |
| --------------------------------------------- | ------- | ------- | ---------- | ---------- | ---------------------------------------------------------------------------------- |
| + [scatter_count](#allOf_i2_scatter_count )   | No      | integer | No         | -          | The amount of scattering that should happen per sample.                            |
| - [merge_distance](#allOf_i2_merge_distance ) | No      | integer | No         | -          | The merge distance for genotype BED files                                          |
| - [dragstr](#allOf_i2_dragstr )               | No      | boolean | No         | -          | Create DragSTR models to be used with HaplotypeCaller                              |
| - [validate](#allOf_i2_validate )             | No      | boolean | No         | -          | Validate the found variants                                                        |
| - [filter](#allOf_i2_filter )                 | No      | boolean | No         | -          | Filter the found variants                                                          |
| - [annotate](#allOf_i2_annotate )             | No      | boolean | No         | -          | Annotate the found variants                                                        |
| - [add_ped](#allOf_i2_add_ped )               | No      | boolean | No         | -          | Add PED INFO header lines to the final VCFs                                        |
| - [gemini](#allOf_i2_gemini )                 | No      | boolean | No         | -          | Create a Gemini databases from the final VCFs                                      |
| - [coverage_fast](#allOf_i2_coverage_fast )   | No      | boolean | No         | -          | Run mosdepth in fast-mode                                                          |
| - [roi](#allOf_i2_roi )                       | No      | string  | No         | -          | Path to the default ROI (regions of interest) BED file to be used for WES analysis |
| - [dbsnp](#allOf_i2_dbsnp )                   | No      | string  | No         | -          | Path to the dbSNP VCF file                                                         |
| - [dbsnp_tbi](#allOf_i2_dbsnp_tbi )           | No      | string  | No         | -          | Path to the index of the dbSNP VCF file                                            |
| - [somalier_sites](#allOf_i2_somalier_sites ) | No      | string  | No         | -          | Path to the VCF file with sites for Somalier to use                                |

### <a name="allOf_i2_scatter_count"></a>3.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > scatter_count`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | Yes       |
| **Default**  | `40`      |

**Description:** The amount of scattering that should happen per sample.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

### <a name="allOf_i2_merge_distance"></a>3.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > merge_distance`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |
| **Default**  | `10000`   |

**Description:** The merge distance for genotype BED files

### <a name="allOf_i2_dragstr"></a>3.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > dragstr`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Create DragSTR models to be used with HaplotypeCaller

### <a name="allOf_i2_validate"></a>3.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > validate`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Validate the found variants

### <a name="allOf_i2_filter"></a>3.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > filter`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Filter the found variants

### <a name="allOf_i2_annotate"></a>3.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > annotate`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Annotate the found variants

### <a name="allOf_i2_add_ped"></a>3.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > add_ped`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Add PED INFO header lines to the final VCFs

### <a name="allOf_i2_gemini"></a>3.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > gemini`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Create a Gemini databases from the final VCFs

### <a name="allOf_i2_coverage_fast"></a>3.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > coverage_fast`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Run mosdepth in fast-mode

### <a name="allOf_i2_roi"></a>3.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > roi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the default ROI (regions of interest) BED file to be used for WES analysis

| Restrictions                      |                                                                                                 |
| --------------------------------- | ----------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.bed(\.gz)?$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.bed%28%5C.gz%29%3F%24) |

### <a name="allOf_i2_dbsnp"></a>3.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > dbsnp`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the dbSNP VCF file

| Restrictions                      |                                                                                     |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.vcf\.gz$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.vcf%5C.gz%24) |

### <a name="allOf_i2_dbsnp_tbi"></a>3.12. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > dbsnp_tbi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the index of the dbSNP VCF file

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.tbi$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.tbi%24) |

### <a name="allOf_i2_somalier_sites"></a>3.13. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Pipeline specific parameters > somalier_sites`

|              |                                                                        |
| ------------ | ---------------------------------------------------------------------- |
| **Type**     | `string`                                                               |
| **Required** | No                                                                     |
| **Format**   | `file-path`                                                            |
| **Default**  | `"https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz"` |

**Description:** Path to the VCF file with sites for Somalier to use

| Restrictions                      |                                                                                 |
| --------------------------------- | ------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.vcf\.gz``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.vcf%5C.gz) |

## <a name="allOf_i3"></a>4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > institutional_config_options`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Defined in**            | #/definitions/institutional_config_options                                |

**Description:** Parameters used to describe centralised config profiles. These should not be edited.

| Property                                                              | Pattern | Type   | Deprecated | Definition | Title/Description                         |
| --------------------------------------------------------------------- | ------- | ------ | ---------- | ---------- | ----------------------------------------- |
| - [custom_config_version](#allOf_i3_custom_config_version )           | No      | string | No         | -          | Git commit id for Institutional configs.  |
| - [custom_config_base](#allOf_i3_custom_config_base )                 | No      | string | No         | -          | Base directory for Institutional configs. |
| - [config_profile_name](#allOf_i3_config_profile_name )               | No      | string | No         | -          | Institutional config name.                |
| - [config_profile_description](#allOf_i3_config_profile_description ) | No      | string | No         | -          | Institutional config description.         |
| - [config_profile_contact](#allOf_i3_config_profile_contact )         | No      | string | No         | -          | Institutional config contact information. |
| - [config_profile_url](#allOf_i3_config_profile_url )                 | No      | string | No         | -          | Institutional config URL link.            |

### <a name="allOf_i3_custom_config_version"></a>4.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > custom_config_version`

|              |            |
| ------------ | ---------- |
| **Type**     | `string`   |
| **Required** | No         |
| **Default**  | `"master"` |

**Description:** Git commit id for Institutional configs.

### <a name="allOf_i3_custom_config_base"></a>4.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > custom_config_base`

|              |                                                              |
| ------------ | ------------------------------------------------------------ |
| **Type**     | `string`                                                     |
| **Required** | No                                                           |
| **Default**  | `"https://raw.githubusercontent.com/nf-core/configs/master"` |

**Description:** Base directory for Institutional configs.

### <a name="allOf_i3_config_profile_name"></a>4.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_name`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Institutional config name.

### <a name="allOf_i3_config_profile_description"></a>4.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_description`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Institutional config description.

### <a name="allOf_i3_config_profile_contact"></a>4.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_contact`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Institutional config contact information.

### <a name="allOf_i3_config_profile_url"></a>4.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Institutional config options > config_profile_url`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Institutional config URL link.

## <a name="allOf_i4"></a>5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > max_job_request_options`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Defined in**            | #/definitions/max_job_request_options                                     |

**Description:** Set the top limit for requested resources for any single job.

| Property                              | Pattern | Type    | Deprecated | Definition | Title/Description                                                  |
| ------------------------------------- | ------- | ------- | ---------- | ---------- | ------------------------------------------------------------------ |
| - [max_cpus](#allOf_i4_max_cpus )     | No      | integer | No         | -          | Maximum number of CPUs that can be requested for any single job.   |
| - [max_memory](#allOf_i4_max_memory ) | No      | string  | No         | -          | Maximum amount of memory that can be requested for any single job. |
| - [max_time](#allOf_i4_max_time )     | No      | string  | No         | -          | Maximum amount of time that can be requested for any single job.   |

### <a name="allOf_i4_max_cpus"></a>5.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Max job request options > max_cpus`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |
| **Default**  | `16`      |

**Description:** Maximum number of CPUs that can be requested for any single job.

### <a name="allOf_i4_max_memory"></a>5.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Max job request options > max_memory`

|              |            |
| ------------ | ---------- |
| **Type**     | `string`   |
| **Required** | No         |
| **Default**  | `"128.GB"` |

**Description:** Maximum amount of memory that can be requested for any single job.

| Restrictions                      |                                                                                                                                                    |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\d+(\.\d+)?\.?\s*(K\|M\|G\|T)?B$``` [Test](https://regex101.com/?regex=%5E%5Cd%2B%28%5C.%5Cd%2B%29%3F%5C.%3F%5Cs%2A%28K%7CM%7CG%7CT%29%3FB%24) |

### <a name="allOf_i4_max_time"></a>5.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Max job request options > max_time`

|              |           |
| ------------ | --------- |
| **Type**     | `string`  |
| **Required** | No        |
| **Default**  | `"240.h"` |

**Description:** Maximum amount of time that can be requested for any single job.

| Restrictions                      |                                                                                                                                            |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^(\d+\.?\s*(s\|m\|h\|day)\s*)+$``` [Test](https://regex101.com/?regex=%5E%28%5Cd%2B%5C.%3F%5Cs%2A%28s%7Cm%7Ch%7Cday%29%5Cs%2A%29%2B%24) |

## <a name="allOf_i5"></a>6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > generic_options`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Defined in**            | #/definitions/generic_options                                             |

**Description:** Less common options for the pipeline, typically set in a config file.

| Property                                                                | Pattern | Type             | Deprecated | Definition | Title/Description                                                                            |
| ----------------------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | -------------------------------------------------------------------------------------------- |
| - [help](#allOf_i5_help )                                               | No      | boolean          | No         | -          | Display help text.                                                                           |
| - [version](#allOf_i5_version )                                         | No      | boolean          | No         | -          | Display version and exit.                                                                    |
| - [publish_dir_mode](#allOf_i5_publish_dir_mode )                       | No      | enum (of string) | No         | -          | Method used to save pipeline results to output directory.                                    |
| - [email_on_fail](#allOf_i5_email_on_fail )                             | No      | string           | No         | -          | Email address for completion summary, only when pipeline fails.                              |
| - [plaintext_email](#allOf_i5_plaintext_email )                         | No      | boolean          | No         | -          | Send plain-text email instead of HTML.                                                       |
| - [max_multiqc_email_size](#allOf_i5_max_multiqc_email_size )           | No      | string           | No         | -          | File size limit when attaching MultiQC reports to summary emails.                            |
| - [monochrome_logs](#allOf_i5_monochrome_logs )                         | No      | boolean          | No         | -          | Do not use coloured log outputs.                                                             |
| - [hook_url](#allOf_i5_hook_url )                                       | No      | string           | No         | -          | Incoming hook URL for messaging service                                                      |
| - [multiqc_title](#allOf_i5_multiqc_title )                             | No      | string           | No         | -          | MultiQC report title. Printed as page header, used for filename if not otherwise specified.  |
| - [multiqc_config](#allOf_i5_multiqc_config )                           | No      | string           | No         | -          | Custom config file to supply to MultiQC.                                                     |
| - [multiqc_logo](#allOf_i5_multiqc_logo )                               | No      | string           | No         | -          | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file |
| - [multiqc_methods_description](#allOf_i5_multiqc_methods_description ) | No      | string           | No         | -          | Custom MultiQC yaml file containing HTML including a methods description.                    |
| - [tracedir](#allOf_i5_tracedir )                                       | No      | string           | No         | -          | Directory to keep pipeline Nextflow logs and reports.                                        |
| - [validate_params](#allOf_i5_validate_params )                         | No      | boolean          | No         | -          | Boolean whether to validate parameters against the schema at runtime                         |
| - [show_hidden_params](#allOf_i5_show_hidden_params )                   | No      | boolean          | No         | -          | Show all params when using \`--help\`                                                        |

### <a name="allOf_i5_help"></a>6.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > help`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Display help text.

### <a name="allOf_i5_version"></a>6.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > version`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Display version and exit.

### <a name="allOf_i5_publish_dir_mode"></a>6.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > publish_dir_mode`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | No                 |
| **Default**  | `"copy"`           |

**Description:** Method used to save pipeline results to output directory.

Must be one of:
* "symlink"
* "rellink"
* "link"
* "copy"
* "copyNoFollow"
* "move"

### <a name="allOf_i5_email_on_fail"></a>6.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > email_on_fail`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Email address for completion summary, only when pipeline fails.

| Restrictions                      |                                                                                                                                                                                                                   |
| --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$``` [Test](https://regex101.com/?regex=%5E%28%5Ba-zA-Z0-9_%5C-%5C.%5D%2B%29%40%28%5Ba-zA-Z0-9_%5C-%5C.%5D%2B%29%5C.%28%5Ba-zA-Z%5D%7B2%2C5%7D%29%24) |

### <a name="allOf_i5_plaintext_email"></a>6.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > plaintext_email`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Send plain-text email instead of HTML.

### <a name="allOf_i5_max_multiqc_email_size"></a>6.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > max_multiqc_email_size`

|              |           |
| ------------ | --------- |
| **Type**     | `string`  |
| **Required** | No        |
| **Default**  | `"25.MB"` |

**Description:** File size limit when attaching MultiQC reports to summary emails.

| Restrictions                      |                                                                                                                                                    |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\d+(\.\d+)?\.?\s*(K\|M\|G\|T)?B$``` [Test](https://regex101.com/?regex=%5E%5Cd%2B%28%5C.%5Cd%2B%29%3F%5C.%3F%5Cs%2A%28K%7CM%7CG%7CT%29%3FB%24) |

### <a name="allOf_i5_monochrome_logs"></a>6.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > monochrome_logs`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Do not use coloured log outputs.

### <a name="allOf_i5_hook_url"></a>6.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > hook_url`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Incoming hook URL for messaging service

### <a name="allOf_i5_multiqc_title"></a>6.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_title`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** MultiQC report title. Printed as page header, used for filename if not otherwise specified.

### <a name="allOf_i5_multiqc_config"></a>6.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_config`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Custom config file to supply to MultiQC.

### <a name="allOf_i5_multiqc_logo"></a>6.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_logo`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file

### <a name="allOf_i5_multiqc_methods_description"></a>6.12. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > multiqc_methods_description`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Custom MultiQC yaml file containing HTML including a methods description.

### <a name="allOf_i5_tracedir"></a>6.13. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > tracedir`

|              |                                    |
| ------------ | ---------------------------------- |
| **Type**     | `string`                           |
| **Required** | No                                 |
| **Default**  | `"${params.outdir}/pipeline_info"` |

**Description:** Directory to keep pipeline Nextflow logs and reports.

### <a name="allOf_i5_validate_params"></a>6.14. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > validate_params`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |
| **Default**  | `true`    |

**Description:** Boolean whether to validate parameters against the schema at runtime

### <a name="allOf_i5_show_hidden_params"></a>6.15. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Generic options > show_hidden_params`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Show all params when using `--help`

## <a name="allOf_i6"></a>7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > annotation_parameters`

|                           |                                                                           |
| ------------------------- | ------------------------------------------------------------------------- |
| **Type**                  | `object`                                                                  |
| **Required**              | No                                                                        |
| **Additional properties** | [[Any type: allowed]](# "Additional Properties of any type are allowed.") |
| **Default**               | `""`                                                                      |
| **Defined in**            | #/definitions/annotation_parameters                                       |

**Description:** Parameters to configure Ensembl VEP and VCFanno

| Property                                              | Pattern | Type    | Deprecated | Definition | Title/Description                                                                                           |
| ----------------------------------------------------- | ------- | ------- | ---------- | ---------- | ----------------------------------------------------------------------------------------------------------- |
| - [species](#allOf_i6_species )                       | No      | string  | No         | -          | The species of the samples                                                                                  |
| - [vep_merged](#allOf_i6_vep_merged )                 | No      | boolean | No         | -          | Specify if the VEP cache is a merged cache                                                                  |
| - [vep_cache](#allOf_i6_vep_cache )                   | No      | string  | No         | -          | The path to the VEP cache                                                                                   |
| - [vep_dbnsfp](#allOf_i6_vep_dbnsfp )                 | No      | boolean | No         | -          | Use the dbNSFP plugin with Ensembl VEP                                                                      |
| - [vep_spliceai](#allOf_i6_vep_spliceai )             | No      | boolean | No         | -          | Use the SpliceAI plugin with Ensembl VEP                                                                    |
| - [vep_spliceregion](#allOf_i6_vep_spliceregion )     | No      | boolean | No         | -          | Use the SpliceRegion plugin with Ensembl VEP                                                                |
| - [vep_mastermind](#allOf_i6_vep_mastermind )         | No      | boolean | No         | -          | Use the Mastermind plugin with Ensembl VEP                                                                  |
| - [vep_maxentscan](#allOf_i6_vep_maxentscan )         | No      | boolean | No         | -          | Use the MaxEntScan plugin with Ensembl VEP                                                                  |
| - [vep_eog](#allOf_i6_vep_eog )                       | No      | boolean | No         | -          | Use the custom EOG annotation with Ensembl VEP                                                              |
| - [vep_version](#allOf_i6_vep_version )               | No      | string  | No         | -          | The version of the VEP tool to be used                                                                      |
| - [vep_cache_version](#allOf_i6_vep_cache_version )   | No      | integer | No         | -          | The version of the VEP cache to be used                                                                     |
| - [dbnsfp](#allOf_i6_dbnsfp )                         | No      | string  | No         | -          | Path to the dbSNFP file                                                                                     |
| - [dbnsfp_tbi](#allOf_i6_dbnsfp_tbi )                 | No      | string  | No         | -          | Path to the index of the dbSNFP file                                                                        |
| - [spliceai_indel](#allOf_i6_spliceai_indel )         | No      | string  | No         | -          | Path to the VCF containing indels for spliceAI                                                              |
| - [spliceai_indel_tbi](#allOf_i6_spliceai_indel_tbi ) | No      | string  | No         | -          | Path to the index of the VCF containing indels for spliceAI                                                 |
| - [spliceai_snv](#allOf_i6_spliceai_snv )             | No      | string  | No         | -          | Path to the VCF containing SNVs for spliceAI                                                                |
| - [spliceai_snv_tbi](#allOf_i6_spliceai_snv_tbi )     | No      | string  | No         | -          | Path to the index of the VCF containing SNVs for spliceAI                                                   |
| - [mastermind](#allOf_i6_mastermind )                 | No      | string  | No         | -          | Path to the VCF for Mastermind                                                                              |
| - [mastermind_tbi](#allOf_i6_mastermind_tbi )         | No      | string  | No         | -          | Path to the index of the VCF for Mastermind                                                                 |
| - [maxentscan](#allOf_i6_maxentscan )                 | No      | string  | No         | -          | Path to the MaxEntScan executable directory                                                                 |
| - [eog](#allOf_i6_eog )                               | No      | string  | No         | -          | Path to the VCF containing EOG annotations                                                                  |
| - [eog_tbi](#allOf_i6_eog_tbi )                       | No      | string  | No         | -          | Path to the index of the VCF containing EOG annotations                                                     |
| - [vcfanno](#allOf_i6_vcfanno )                       | No      | boolean | No         | -          | Run annotations with vcfanno                                                                                |
| - [vcfanno_config](#allOf_i6_vcfanno_config )         | No      | string  | No         | -          | The path to the VCFanno config TOML                                                                         |
| - [vcfanno_lua](#allOf_i6_vcfanno_lua )               | No      | string  | No         | -          | The path to a Lua script to be used in VCFanno                                                              |
| - [vcfanno_resources](#allOf_i6_vcfanno_resources )   | No      | string  | No         | -          | A comma-seperated list of resource files for VCFanno, please also supply their indices using this parameter |

### <a name="allOf_i6_species"></a>7.1. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > species`

|              |                  |
| ------------ | ---------------- |
| **Type**     | `string`         |
| **Required** | No               |
| **Default**  | `"homo_sapiens"` |

**Description:** The species of the samples

| Restrictions                      |                                                                         |
| --------------------------------- | ----------------------------------------------------------------------- |
| **Must match regular expression** | ```^[a-z_]*$``` [Test](https://regex101.com/?regex=%5E%5Ba-z_%5D%2A%24) |

### <a name="allOf_i6_vep_merged"></a>7.2. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_merged`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |
| **Default**  | `true`    |

**Description:** Specify if the VEP cache is a merged cache

### <a name="allOf_i6_vep_cache"></a>7.3. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_cache`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |
| **Format**   | `path`   |
| **Default**  | `"None"` |

**Description:** The path to the VEP cache

### <a name="allOf_i6_vep_dbnsfp"></a>7.4. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_dbnsfp`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Use the dbNSFP plugin with Ensembl VEP

### <a name="allOf_i6_vep_spliceai"></a>7.5. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_spliceai`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Use the SpliceAI plugin with Ensembl VEP

### <a name="allOf_i6_vep_spliceregion"></a>7.6. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_spliceregion`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Use the SpliceRegion plugin with Ensembl VEP

### <a name="allOf_i6_vep_mastermind"></a>7.7. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_mastermind`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Use the Mastermind plugin with Ensembl VEP

### <a name="allOf_i6_vep_maxentscan"></a>7.8. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_maxentscan`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Use the MaxEntScan plugin with Ensembl VEP

### <a name="allOf_i6_vep_eog"></a>7.9. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_eog`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Use the custom EOG annotation with Ensembl VEP

### <a name="allOf_i6_vep_version"></a>7.10. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_version`

|              |           |
| ------------ | --------- |
| **Type**     | `string`  |
| **Required** | No        |
| **Default**  | `"105.0"` |

**Description:** The version of the VEP tool to be used

| Restrictions                      |                                                                                   |
| --------------------------------- | --------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\d\d\d?\.\d``` [Test](https://regex101.com/?regex=%5E%5Cd%5Cd%5Cd%3F%5C.%5Cd) |

### <a name="allOf_i6_vep_cache_version"></a>7.11. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vep_cache_version`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |
| **Default**  | `105`     |

**Description:** The version of the VEP cache to be used

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

### <a name="allOf_i6_dbnsfp"></a>7.12. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > dbnsfp`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the dbSNFP file

| Restrictions                      |                                                                         |
| --------------------------------- | ----------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.gz$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.gz%24) |

### <a name="allOf_i6_dbnsfp_tbi"></a>7.13. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > dbnsfp_tbi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the index of the dbSNFP file

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.csi$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.csi%24) |

### <a name="allOf_i6_spliceai_indel"></a>7.14. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_indel`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the VCF containing indels for spliceAI

| Restrictions                      |                                                                                     |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.vcf\.gz$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.vcf%5C.gz%24) |

### <a name="allOf_i6_spliceai_indel_tbi"></a>7.15. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_indel_tbi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the index of the VCF containing indels for spliceAI

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.tbi$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.tbi%24) |

### <a name="allOf_i6_spliceai_snv"></a>7.16. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_snv`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the VCF containing SNVs for spliceAI

| Restrictions                      |                                                                                     |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.vcf\.gz$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.vcf%5C.gz%24) |

### <a name="allOf_i6_spliceai_snv_tbi"></a>7.17. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > spliceai_snv_tbi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the index of the VCF containing SNVs for spliceAI

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.tbi$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.tbi%24) |

### <a name="allOf_i6_mastermind"></a>7.18. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > mastermind`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the VCF for Mastermind

| Restrictions                      |                                                                                     |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.vcf\.gz$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.vcf%5C.gz%24) |

### <a name="allOf_i6_mastermind_tbi"></a>7.19. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > mastermind_tbi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the index of the VCF for Mastermind

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.tbi$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.tbi%24) |

### <a name="allOf_i6_maxentscan"></a>7.20. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > maxentscan`

|              |                  |
| ------------ | ---------------- |
| **Type**     | `string`         |
| **Required** | No               |
| **Format**   | `directory-path` |

**Description:** Path to the MaxEntScan executable directory

### <a name="allOf_i6_eog"></a>7.21. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > eog`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the VCF containing EOG annotations

| Restrictions                      |                                                                                     |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.vcf\.gz$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.vcf%5C.gz%24) |

### <a name="allOf_i6_eog_tbi"></a>7.22. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > eog_tbi`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** Path to the index of the VCF containing EOG annotations

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.tbi$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.tbi%24) |

### <a name="allOf_i6_vcfanno"></a>7.23. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |

**Description:** Run annotations with vcfanno

### <a name="allOf_i6_vcfanno_config"></a>7.24. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno_config`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** The path to the VCFanno config TOML

| Restrictions                      |                                                                             |
| --------------------------------- | --------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.toml$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.toml%24) |

### <a name="allOf_i6_vcfanno_lua"></a>7.25. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno_lua`

|              |             |
| ------------ | ----------- |
| **Type**     | `string`    |
| **Required** | No          |
| **Format**   | `file-path` |

**Description:** The path to a Lua script to be used in VCFanno

| Restrictions                      |                                                                           |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Must match regular expression** | ```^\S+\.lua$``` [Test](https://regex101.com/?regex=%5E%5CS%2B%5C.lua%24) |

### <a name="allOf_i6_vcfanno_resources"></a>7.26. Property `CenterForMedicalGeneticsGhent/nf-cmgg-germline pipeline parameters > allOf > Annotation parameters > vcfanno_resources`

|              |                  |
| ------------ | ---------------- |
| **Type**     | `string`         |
| **Required** | No               |
| **Format**   | `directory-path` |

**Description:** A comma-seperated list of resource files for VCFanno, please also supply their indices using this parameter

----------------------------------------------------------------------------------------------------------------------------
Generated using [json-schema-for-humans](https://github.com/coveooss/json-schema-for-humans) on 2023-04-25 at 10:03:53 +0200
