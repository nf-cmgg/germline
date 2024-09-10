#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/germline
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getGenomeAttribute } from './subworkflows/local/utils_cmgg_germline_pipeline'

params.fasta                = getGenomeAttribute('fasta', params.genomes, params.genome)
params.fai                  = getGenomeAttribute('fai', params.genomes, params.genome)
params.dict                 = getGenomeAttribute('dict', params.genomes, params.genome)
params.strtablefile         = getGenomeAttribute('strtablefile', params.genomes, params.genome)
params.sdf                  = getGenomeAttribute('sdf', params.genomes, params.genome)
params.dbsnp                = getGenomeAttribute('dbsnp', params.genomes, params.genome)
params.dbsnp_tbi            = getGenomeAttribute('dbsnp_tbi', params.genomes, params.genome)
params.vep_cache            = getGenomeAttribute('vep_cache', params.genomes, params.genome)
params.dbnsfp               = getGenomeAttribute('dbnsfp', params.genomes, params.genome)
params.dbnsfp_tbi           = getGenomeAttribute('dbnsfp_tbi', params.genomes, params.genome)
params.spliceai_indel       = getGenomeAttribute('spliceai_indel', params.genomes, params.genome)
params.spliceai_indel_tbi   = getGenomeAttribute('spliceai_indel_tbi', params.genomes, params.genome)
params.spliceai_snv         = getGenomeAttribute('spliceai_snv', params.genomes, params.genome)
params.spliceai_snv_tbi     = getGenomeAttribute('spliceai_snv_tbi', params.genomes, params.genome)
params.mastermind           = getGenomeAttribute('mastermind', params.genomes, params.genome)
params.mastermind_tbi       = getGenomeAttribute('mastermind_tbi', params.genomes, params.genome)
params.eog                  = getGenomeAttribute('eog', params.genomes, params.genome)
params.eog_tbi              = getGenomeAttribute('eog_tbi', params.genomes, params.genome)
params.alphamissense        = getGenomeAttribute('alphamissense', params.genomes, params.genome)
params.alphamissense_tbi    = getGenomeAttribute('alphamissense_tbi', params.genomes, params.genome)
params.vcfanno_resources    = getGenomeAttribute('vcfanno_resources', params.genomes, params.genome)
params.vcfanno_config       = getGenomeAttribute('vcfanno_config', params.genomes, params.genome)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GERMLINE                } from './workflows/germline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_cmgg_germline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_cmgg_germline_pipeline'

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//

workflow NFCMGG_GERMLINE {

    take:
    samplesheet     // channel: samplesheet read in from --input
    pipeline_params // the parameters used for this pipeline
    multiqc_logo    // string: the path to the multiqc logo

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GERMLINE (
        // Input channels
        samplesheet,

        // File inputs
        pipeline_params.fasta,
        pipeline_params.fai,
        pipeline_params.dict,
        pipeline_params.strtablefile,
        pipeline_params.sdf,
        pipeline_params.dbsnp,
        pipeline_params.dbsnp_tbi,
        pipeline_params.vep_cache,
        pipeline_params.dbnsfp,
        pipeline_params.dbnsfp_tbi,
        pipeline_params.spliceai_indel,
        pipeline_params.spliceai_indel_tbi,
        pipeline_params.spliceai_snv,
        pipeline_params.spliceai_snv_tbi,
        pipeline_params.mastermind,
        pipeline_params.mastermind_tbi,
        pipeline_params.eog,
        pipeline_params.eog_tbi,
        pipeline_params.alphamissense,
        pipeline_params.alphamissense_tbi,
        pipeline_params.vcfanno_resources,
        pipeline_params.vcfanno_config,
        pipeline_params.multiqc_config,
        multiqc_logo,
        pipeline_params.multiqc_methods_description,
        pipeline_params.roi,
        pipeline_params.somalier_sites,
        pipeline_params.vcfanno_lua,
        pipeline_params.updio_common_cnvs,
        pipeline_params.automap_repeats,
        pipeline_params.automap_panel,
        pipeline_params.outdir,
        GlobalVariables.pedFiles,

        // Boolean inputs
        pipeline_params.dragstr,
        pipeline_params.annotate,
        pipeline_params.vcfanno,
        pipeline_params.only_call,
        pipeline_params.only_merge,
        pipeline_params.filter,
        pipeline_params.normalize,
        pipeline_params.add_ped,
        pipeline_params.gemini,
        pipeline_params.validate,
        pipeline_params.updio,
        pipeline_params.automap,
        pipeline_params.vep_dbnsfp,
        pipeline_params.vep_spliceai,
        pipeline_params.vep_mastermind,
        pipeline_params.vep_eog,
        pipeline_params.vep_alphamissense,

        // Value inputs
        pipeline_params.genome,
        pipeline_params.species,
        pipeline_params.vep_cache_version,
        pipeline_params.vep_chunk_size,
        pipeline_params.scatter_count,
        pipeline_params.callers.tokenize(",")
    )

    emit:
    multiqc_report = GERMLINE.out.multiqc_report // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

    if(params.dbsnp_tbi && !params.dbsnp){
        error("Please specify the dbsnp VCF with --dbsnp VCF")
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

    callers = params.callers.tokenize(",")
    callers.each { caller ->
        if(!(caller in GlobalVariables.availableCallers)) { error("\"${caller}\" is not a supported callers please use one or more of these instead: ${GlobalVariables.availableCallers}")}
    }

    if (params.output_suffix && callers.size() > 1) {
        error("Cannot use --output_suffix with more than one caller")
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def multiqc_logo = params.multiqc_logo   ?: "$projectDir/assets/CMGG_logo.png"

    print(params.genomes)
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.ped,
        params.genomes,
        params.genome
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCMGG_GERMLINE (
        PIPELINE_INITIALISATION.out.samplesheet,
        params,
        multiqc_logo
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
        params.hook_url,
        NFCMGG_GERMLINE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
