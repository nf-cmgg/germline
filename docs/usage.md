# nf-cmgg/germline: Usage

> _Documentation of pipeline parameters can be found in the [parameters documentation](./parameters.md)_

## Samplesheet input

You will need to create a samplesheet with information with the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It can be either a CSV, TSV, JSON or YAML file.

```bash
--input '[path to samplesheet file]'
```

### Watch for files in a directory

When the `--watchdir` parameter has been given, the pipeline will automatically check for all files in the samplesheet that have the `watch:` prefix in the given directory. An example for watching CRAM files:

```csv title="samplesheet.csv"
sample,cram,crai
SAMPLE_1,watch:INPUT.cram,watch:INPUT.cram.crai
```

The files `INPUT.cram` and `INPUT.cram.crai` will now be watched for recursively in the watch directory.

### Example of the samplesheet

Below is an example of how the samplesheet could look like in the three formats.

<!-- prettier-ignore -->
!!!note
    The order and presence of the fields is not set, you can arrange/remove these as you see fit. The only required fields are `sample` and `cram`.

#### CSV

```csv title="samplesheet.csv"
sample,family,cram,crai
SAMPLE_1,FAMILY_1,SAMPLE_1.cram,SAMPLE_1.crai
SAMPLE_2,FAMILY_1,SAMPLE_2.cram,SAMPLE_2.crai
SAMPLE_3,,SAMPLE_3.cram,
```

#### TSV

```tsv title="samplesheet.tsv"
sample    family    cram   crai
SAMPLE_1  FAMILY_1  SAMPLE_1.cram SAMPLE_1.crai
SAMPLE_2  FAMILY_1  SAMPLE_2.cram SAMPLE_2.crai
SAMPLE_3    SAMPLE_3.cram
```

#### YAML/YML

```yaml title="samplesheet.yaml"
- sample: SAMPLE_1
  family: FAMILY_1
  cram: SAMPLE_1.cram
  crai: SAMPLE_1.crai
- sample: SAMPLE_2
  family: FAMILY_1
  cram: SAMPLE_2.cram
  crai: SAMPLE_2.crai
- sample: SAMPLE_3
  cram: SAMPLE_3.cram
```

#### JSON

```json title="samplesheet.json"
[
  {
    "sample": "SAMPLE_1",
    "family": "FAMILY_1",
    "cram": "SAMPLE_1.cram",
    "crai": "SAMPLE_1.crai"
  },
  {
    "sample": "SAMPLE_2",
    "family": "FAMILY_1",
    "cram": "SAMPLE_2.cram",
    "crai": "SAMPLE_2.crai"
  },
  {
    "sample": "SAMPLE_3",
    "cram": "SAMPLE_3.cram"
  }
]
```

### Full samplesheet

The samplesheet can have following columns:

| Column           | Description                                                                                                                                                                                                                                                                                                                                                           |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`         | MANDATORY - Custom sample name. This entry has to be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`).                                                                                                                                                                  |
| `family`         | OPTIONAL - The family ID of the specified sample. This field is optional, as the family id can also be extracted from the `ped` file. If no `ped` file and `family` ID are supplied, the `family` ID defaults to the `sample` ID (which means that the resulting VCF will be single-sample). Spaces in family names are automatically converted to underscores (`_`). |
| `cram`           | MANDATORY - Full path to CRAM file to call variants from. File has to have the extension `.cram`                                                                                                                                                                                                                                                                      |
| `crai`           | OPTIONAL - Full path to CRAM index file. File has to have the extension `.crai`.                                                                                                                                                                                                                                                                                      |
| `ped`            | OPTIONAL - Full path to PED file containing the relational information between samples in the same family. File has to have the extension `.ped`.                                                                                                                                                                                                                     |
| `truth_vcf`      | OPTIONAL - Full path to the VCF containing all the truth variants of the current sample. The validation subworkflow will be run when this file is supplied and the `--validate true` flag has been given. File has to have the extension `.vcf.gz`                                                                                                                    |
| `truth_tbi`      | OPTIONAL - Full path to the index of the truth VCF. This file can either be supplied by the user or generated by the pipeline. File has to have the extensions `.tbi`                                                                                                                                                                                                 |
| `truth_bed`      | OPTIONAL - Full path to the BED file containing the golden truth regions in the `truth_vcf` file. File has to have the extensions `.bed`                                                                                                                                                                                                                              |
| `roi`            | OPTIONAL - Full path to a BED file containing the regions of interest for the current sample to call on. When this file is given, the pipeline will run this sample in WES mode. (The flag `--roi <path>` can also be given to run WES mode for all samples using the file specified by the flag) File has to have the extension `.bed` or `.bed.gz`.                 |
| `vardict_min_af` | OPTIONAL - The minimum AF value to use for the vardict variant caller (`--callers vardict`). This can be set in the samplesheet when it differs for all samples. A default can be set using the `--vardict_min_af` parameter (whichs defaults to 0.1)                                                                                                                 |

<!-- prettier-ignore -->
!!!note
    The `sample` fields has to contain the same value when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. Either the `ped` or `family` field can be used to specify the family name. The pipeline automatically extracts the family id from the `ped` file if the `family` field is empty. The `family` is used to specify on which samples the joint-genotyping should be performed. If neither the `ped` or `family` fields are used, the pipeline will default to a single-sample family with the sample name as its ID.

This is an example of a working samplesheet used to test this pipeline:

```csv title="samplesheet.csv"
--8<-- "assets/samplesheet.csv"
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-cmgg/germline --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work          #(1)!
results       #(2)!
.nextflow.log #(3)!
...           #(4)!
```

1. Directory containing the nextflow working files

2. Finished results in specified location (defined with --outdir). See [output](./output.md) documentation for more on this.

3. Log file from Nextflow

4. Other nextflow hidden files, eg. history of pipeline runs and old logs.

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

<!-- prettier-ignore -->
!!!warning
    Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-cmgg/germline -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh38'
<...>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline. You can also add the `-latest` argument to your run command to automatically fetch the latest version on every run:

```bash
nextflow pull nf-cmgg/germline -r <version>
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-cmgg/germline releases page](https://github.com/nf-cmgg/germline/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

<!-- prettier-ignore -->
!!!tip
    If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

<!-- prettier-ignore -->
!!!note
    These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.
Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

<!-- prettier-ignore -->
!!!info
    We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  > A profile with a complete configuration for automated testing
  > Includes links to test data so needs no other parameters
- `nf_test`
  > The profile setting the default values for `nf-test`. When running `nf-test` this profile is automatically used.
- `docker`
  > A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  > A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  > A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  > A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  > A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  > A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  > A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file. See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-cmgg/germline/blob/b637c64c2e1eeb1527d481a377f60950c9a114b8/conf/base.config#L17) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
