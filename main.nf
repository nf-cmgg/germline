#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/smallvariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/smallvariants
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getGenomeAttribute } from './subworkflows/local/utils_cmgg_smallvariants_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SMALLVARIANTS           } from './workflows/smallvariants'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_cmgg_smallvariants_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_cmgg_smallvariants_pipeline'
include { getWorkflowVersion      } from './subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated file containing information about the samples in the experiment.
    input: Path

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String

    // Email address for completion summary.
    email: String?

    // Path to a pedigree file for all samples in the run. All relational data will be fetched from this file.
    ped: Path?

    // Path to FASTA genome file.
    fasta: Path = getGenomeAttribute('fasta', params.genomes, params.genome)

    // Path to FASTA genome index file.
    fai: Path? = getGenomeAttribute('fai', params.genomes, params.genome)

    // Path to the sequence dictionary generated from the FASTA reference. This is only used when `haplotypecaller` is one of the specified callers.
    dict: Path? = getGenomeAttribute('dict', params.genomes, params.genome)

    // Path to the STR table file generated from the FASTA reference. This is only used when `--dragstr` has been given.
    strtablefile: Path? = getGenomeAttribute('strtablefile', params.genomes, params.genome)

    // Path to the SDF folder generated from the reference FASTA file. This is only required when using `--validate`.
    sdf: Path? = getGenomeAttribute('sdf', params.genomes, params.genome)

    // Path to the ELFASTA genome file. This is used when `elprep` is part of the callers and will be automatically generated when missing.
    elfasta: Path? = getGenomeAttribute('elfasta', params.genomes, params.genome)

    // Path to the elsites file. This is used when `elprep` is part of the callers.
    elsites: Path?

    // Object for genomes
    genomes: Map = [:]

    // Directory / URL base for iGenomes references.
    igenomes_base: String = '/references'

    // Do not load the iGenomes reference config.
    igenomes_ignore: Boolean = false

    // Path to the MSI baseline VCF file.
    msi_baseline: Path?

    // The amount of scattering that should happen per sample.
    scatter_count: Integer = 40

    // The merge distance for family BED files
    merge_distance: Integer = 100000

    // Create DragSTR models to be used with HaplotypeCaller
    dragstr: Boolean = false

    // Validate the found variants
    validate: Boolean = false

    // Filter the found variants.
    filter: Boolean = false

    // Annotate the found variants using Ensembl VEP.
    annotate: Boolean = false

    // Add PED INFO header lines to the final VCFs.
    add_ped: Boolean = false

    // Create a Gemini databases from the final VCFs.
    gemini: Boolean = false

    // Don't run mosdepth in fast-mode
    mosdepth_slow: Boolean = false

    // Path to the default ROI (regions of interest) BED file to be used for WES analysis.
    roi: Path?

    // Path to the dbSNP VCF file. This will be used to set the variant IDs.
    dbsnp: Path? = getGenomeAttribute('dbsnp', params.genomes, params.genome)

    // Path to the index of the dbSNP VCF file.
    dbsnp_tbi: Path? = getGenomeAttribute('dbsnp_tbi', params.genomes, params.genome)

    // Path to the VCF file with sites for Somalier to use.
    somalier_sites: Path? = getGenomeAttribute('somalier_sites', params.genomes, params.genome)

    // Only call the variants without doing any post-processing.
    only_call: Boolean = false

    // Only run the pipeline until the creation of the genomicsdbs and output them.
    only_merge: Boolean = false

    // Output the genomicsDB together with the joint-genotyped VCF.
    output_genomicsdb: Boolean = false

    // A comma delimited string of the available callers. Current options are: `haplotypecaller` and `vardict`.
    callers: String = 'haplotypecaller'

    // The minimum allele frequency for VarDict when no `vardict_min_af` is supplied in the samplesheet.
    vardict_min_af: Float = 0.1

    // Normalize the variant in the final VCFs.
    normalize: Boolean = false

    // Filter out all variants that don't have the PASS filter for vardict. This only works when `--filter` is also given.
    only_pass: Boolean = false

    // Keep all aditional contigs for calling instead of filtering them out before.
    keep_alt_contigs: Boolean = false

    // Disable UPDio analysis on the final VCFs.
    disable_updio: Boolean = false

    // A TSV file containing common CNVs to be used by UPDio.
    updio_common_cnvs: Path?

    // Disable AutoMap analysis on the final VCFs.
    disable_automap: Boolean = false

    // BED file with repeat regions in the genome.
    automap_repeats: Path?

    // TXT file with gene panel regions to be used by AutoMap.
    automap_panel: Path?

    // BED file with regions to be used by UPDio.
    updio_regions: Path?

    // The panel name of the panel given with --automap_panel.
    automap_panel_name: String = 'cmgg_bio'

    // Perform phasing with HaplotypeCaller.
    hc_phasing: Boolean = false

    // The lowest callable coverage to determine callable regions.
    min_callable_coverage: Integer = 5

    // Don't change this value
    unique_out: String = "v${workflow.manifest.version.replace('.', '_')}_${new java.util.Date().format('yyyy_MM_dd')}"

    // Disable the sequence dictionary validation in HaplotypeCaller
    disable_hc_dict_validation: Boolean = false

    // Don't output the merged CRAM files.
    skip_merged_cram_output: Boolean = false

    // Display version and exit.
    version: Boolean = false

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean = false

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean = false

    // Incoming hook URL for messaging service
    hook_url: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String = 'https://raw.githubusercontent.com/nf-core/test-datasets/'

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String

    // Display the help message.
    help = false

    // Display the full detailed help message.
    help_full: Boolean = false

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean = false

    // The amount of sites per split VCF as input to VEP.
    vep_chunk_size: Integer = 50000

    // The species of the samples.
    species: String = 'homo_sapiens'

    // Specify if the VEP cache is a merged cache.
    vep_merged: Boolean = true

    // The path to the VEP cache.
    vep_cache: Path? = getGenomeAttribute('vep_cache', params.genomes, params.genome)

    // Use the dbNSFP plugin with Ensembl VEP.
    vep_dbnsfp: Boolean = false

    // Use the SpliceAI plugin with Ensembl VEP.
    vep_spliceai: Boolean = false

    // Use the SpliceRegion plugin with Ensembl VEP.
    vep_spliceregion: Boolean = false

    // Use the Mastermind plugin with Ensembl VEP.
    vep_mastermind: Boolean = false

    // Use the MaxEntScan plugin with Ensembl VEP.
    vep_maxentscan: Boolean = false

    // Use the custom EOG annotation with Ensembl VEP.
    vep_eog: Boolean = false

    // Use the AlphaMissense plugin with Ensembl VEP.
    vep_alphamissense: Boolean = false

    // The version of the VEP tool to be used.
    vep_version: Float = 105.0

    // The version of the VEP cache to be used.
    vep_cache_version: Integer = 105

    // Path to the dbSNFP file.
    dbnsfp: Path? = getGenomeAttribute('dbnsfp', params.genomes, params.genome)

    // Path to the index of the dbSNFP file.
    dbnsfp_tbi: Path? = getGenomeAttribute('dbnsfp_tbi', params.genomes, params.genome)

    // Path to the VCF containing indels for spliceAI.
    spliceai_indel: Path? = getGenomeAttribute('spliceai_indel', params.genomes, params.genome)

    // Path to the index of the VCF containing indels for spliceAI.
    spliceai_indel_tbi: Path? = getGenomeAttribute('spliceai_indel_tbi', params.genomes, params.genome)

    // Path to the VCF containing SNVs for spliceAI.
    spliceai_snv: Path? = getGenomeAttribute('spliceai_snv', params.genomes, params.genome)

    // Path to the index of the VCF containing SNVs for spliceAI.
    spliceai_snv_tbi: Path? = getGenomeAttribute('spliceai_snv_tbi', params.genomes, params.genome)

    // Path to the VCF for Mastermind.
    mastermind: Path? = getGenomeAttribute('mastermind', params.genomes, params.genome)

    // Path to the index of the VCF for Mastermind.
    mastermind_tbi: Path? = getGenomeAttribute('mastermind_tbi', params.genomes, params.genome)

    // Path to the TSV for AlphaMissense.
    alphamissense: Path? = getGenomeAttribute('alphamissense', params.genomes, params.genome)

    // Path to the index of the TSV for AlphaMissense.
    alphamissense_tbi: Path? = getGenomeAttribute('alphamissense_tbi', params.genomes, params.genome)

    // Path to the VCF containing EOG annotations.
    eog: Path? = getGenomeAttribute('eog', params.genomes, params.genome)

    // Path to the index of the VCF containing EOG annotations.
    eog_tbi: Path? = getGenomeAttribute('eog_tbi', params.genomes, params.genome)

    // Path to the directory containing the MaxEntScan reference annotations.
    maxentscan: Path? = getGenomeAttribute('maxentscan', params.genomes, params.genome)

    // Run annotations with vcfanno.
    vcfanno: Boolean = false

    // The path to the VCFanno config TOML.
    vcfanno_config: Path? = getGenomeAttribute('vcfanno_config', params.genomes, params.genome)

    // The path to a Lua script to be used in VCFanno.
    vcfanno_lua: Path?

    // A semicolon-seperated list of resource files for VCFanno, please also supply their indices using this parameter.
    vcfanno_resources: String? = getGenomeAttribute('vcfanno_resources', params.genomes, params.genome)

    // Path to the panel of normals VCF file.
    panel_of_normals: Path? = getGenomeAttribute('panel_of_normals', params.genomes, params.genome)

    // Path to the index of the panel of normals VCF file.
    panel_of_normals_tbi: Path? = getGenomeAttribute('panel_of_normals_tbi', params.genomes, params.genome)
}

workflow {

    main:
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // Check for dependencies between parameters
    //

    def List<String> available_callers = ["haplotypecaller", "vardict", "elprep", "mutect2"]

    if(params.dbsnp_tbi && !params.dbsnp){
        error("Please specify the dbsnp VCF with --dbsnp VCF")
    }

    if(params.panel_of_normals_tbi && !params.panel_of_normals){
        error("Please specify the panel of normals VCF with --panel_of_normals VCF")
    }

    if (params.annotate) {
        // Check if a genome is given
        if (!params.genome) { error("A genome should be supplied for annotation (use --genome)") }

        // Check if the VEP versions were given
        if (!params.vep_version) { error("A VEP version should be supplied for annotation (use --vep_version)") }
        if (!params.vep_cache_version) { error("A VEP cache version should be supplied for annotation (use --vep_cache_version)") }

        // Check if a species is entered
        if (!params.species) { error("A species should be supplied for annotation (use --species)") }

        // Check if all vcfanno files are supplied when vcfanno should be used
        if (params.vcfanno && (!params.vcfanno_config || !params.vcfanno_resources)) {
            error("A TOML file and resource files should be supplied when using vcfanno (use --vcfanno_config and --vcfanno_resources)")
        }
    }

    def callers = params.callers.tokenize(",")
    callers.each { caller ->
        if(!(caller in available_callers)) { error("\"${caller}\" is not a supported caller please use one or more of these instead: ${available_callers.join(', ')}") }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def multiqc_logo = params.multiqc_logo ?: file("$projectDir/assets/CMGG_logo.png")

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.ped,
        params.help,
        params.help_full,
        params.show_hidden,
        params.unique_out,
    )

    //
    // WORKFLOW: Run main workflow
    //

    SMALLVARIANTS (
        // Input channels
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.ped_files,

        // File inputs
        params.fasta,
        params.fai,
        params.dict,
        params.elfasta,
        params.strtablefile,
        params.sdf,
        params.dbsnp,
        params.dbsnp_tbi,
        params.vep_cache,
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.spliceai_indel,
        params.spliceai_indel_tbi,
        params.spliceai_snv,
        params.spliceai_snv_tbi,
        params.mastermind,
        params.mastermind_tbi,
        params.eog,
        params.eog_tbi,
        params.alphamissense,
        params.alphamissense_tbi,
        params.maxentscan,
        params.vcfanno_resources,
        params.vcfanno_config,
        params.multiqc_config,
        multiqc_logo,
        params.multiqc_methods_description,
        params.roi,
        params.somalier_sites,
        params.vcfanno_lua,
        params.updio_common_cnvs,
        params.automap_repeats,
        params.automap_panel,
        params.outdir,
        params.elsites,
        params.msi_baseline,
        params.updio_regions,
        params.panel_of_normals,
        params.panel_of_normals_tbi,

        // Boolean inputs
        params.dragstr,
        params.annotate,
        params.vcfanno,
        params.only_call,
        params.only_merge,
        params.filter,
        params.normalize,
        params.add_ped,
        params.gemini,
        params.validate,
        params.disable_updio,
        params.disable_automap,
        params.vep_dbnsfp,
        params.vep_spliceai,
        params.vep_mastermind,
        params.vep_eog,
        params.vep_alphamissense,
        params.vep_maxentscan,
        params.mosdepth_slow,

        // Value inputs
        params.genome,
        params.species,
        params.vep_cache_version,
        params.vep_chunk_size,
        params.scatter_count,
        params.callers.tokenize(",")
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        SMALLVARIANTS.out.multiqc_report.toList()
    )

    publish:
    merged_crams        = SMALLVARIANTS.out.merged_crams
    mosdepth_reports    = SMALLVARIANTS.out.mosdepth_reports
    gvcfs               = SMALLVARIANTS.out.gvcfs.filter { _meta, gvcf, _tbi -> gvcf.startsWith(workflow.workDir) } // Filtering out input GVCFs from the output publishing fixes an issue in the current implementation of the workflow output definitions: https://github.com/nextflow-io/nextflow/issues/5480
    msi                 = SMALLVARIANTS.out.msi
    single_beds         = SMALLVARIANTS.out.single_beds
    perbase_beds        = SMALLVARIANTS.out.perbase_beds
    validation          = SMALLVARIANTS.out.validation
    gvcf_reports        = SMALLVARIANTS.out.gvcf_reports
    genomicsdb          = SMALLVARIANTS.out.genomicsdb
    vcfs                = SMALLVARIANTS.out.vcfs.filter { _meta, vcf, _tbi -> vcf.startsWith(workflow.workDir) } // Filtering out input VCFs from the output publishing fixes an issue in the current implementation of the workflow output definitions: https://github.com/nextflow-io/nextflow/issues/5480
    gemini              = SMALLVARIANTS.out.gemini
    peds                = SMALLVARIANTS.out.peds
    joint_beds          = SMALLVARIANTS.out.joint_beds
    final_reports       = SMALLVARIANTS.out.final_reports
    automap             = SMALLVARIANTS.out.automap.map { meta, dir -> [ meta, files("${dir.toUriString()}/**") ] }.transpose(by:1)
    updio               = SMALLVARIANTS.out.updio
    multiqc             = SMALLVARIANTS.out.multiqc_report
    multiqc_data        = SMALLVARIANTS.out.multiqc_data
}

output {
    merged_crams {
        path { meta, cram, crai ->
            cram >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.cram"
            crai >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.cram.crai"
        }
        enabled !params.skip_merged_cram_output
    }
    mosdepth_reports { path { meta, report ->
        report >> "${meta.family}/${meta.id}_${params.unique_out}/mosdepth/${report.name}"
    } }
    gvcfs { path { meta, gvcf, tbi ->
        gvcf >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.${meta.caller}.g.vcf.gz"
        tbi >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.${meta.caller}.g.vcf.gz.tbi"
    } }
    msi { path { meta, msi ->
        msi >> "${meta.family}/${meta.id}_${params.unique_out}/msi/${msi.name}"
    } }
    single_beds { path { meta, bed ->
        bed >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.bed"
    } }
    perbase_beds { path { meta, bed, csi ->
        bed >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.per-base.bed.gz"
        csi >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.per-base.bed.gz.csi"
    } }
    validation { path { meta, report ->
        report >> "${meta.family}/${meta.id}_${params.unique_out}/validation/${meta.caller}/${report.name}"
    } }
    gvcf_reports { path { meta, report ->
        report >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.${meta.caller}.bcftools_stats.txt"
    }}
    genomicsdb {
        enabled (params.output_genomicsdb || params.only_merge)
        path { meta, genomicsdb ->
            genomicsdb >> "${meta.family}/output_${params.unique_out}/${meta.id}_${meta.caller}_genomicsdb"
        }
    }
    vcfs { path { meta, vcf, tbi ->
        vcf >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.vcf.gz"
        tbi >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.vcf.gz.tbi"
    } }
    gemini { path { meta, db ->
        db >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.db"
    } }
    peds { path { meta, ped ->
        ped >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.ped"
    } }
    joint_beds { path { meta, bed ->
        bed >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.bed"
    } }
    final_reports { path { meta, report ->
        report >> "${meta.family}/qc_${params.unique_out}/${report.name}"
    } }
    automap { path { meta, automap_file ->
        automap_file >> "${meta.family}/output_${params.unique_out}/automap/${meta.caller}/${automap_file.name}"
    } }
    updio { path { meta, updio ->
        updio >> "${meta.family}/output_${params.unique_out}/updio/${meta.caller}${params.updio_regions ? '.filtered' : ''}"
    } }
    multiqc { path { _meta, report ->
        report >> "${params.unique_out}/multiqc_report.html"
    } }
    multiqc_data { path { _meta, _folder ->
        "${params.unique_out}/"
    } }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
