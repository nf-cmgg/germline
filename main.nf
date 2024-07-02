#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/germline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getGenomeAttribute      } from './subworkflows/local/utils_cmgg_germline_pipeline'

params.fasta                = getGenomeAttribute('fasta')
params.fai                  = getGenomeAttribute('fai')
params.dict                 = getGenomeAttribute('dict')
params.strtablefile         = getGenomeAttribute('strtablefile')
params.sdf                  = getGenomeAttribute('sdf')
params.dbsnp                = getGenomeAttribute('dbsnp')
params.dbsnp_tbi            = getGenomeAttribute('dbsnp_tbi')
params.vep_cache            = getGenomeAttribute('vep_cache')
params.dbnsfp               = getGenomeAttribute('dbnsfp')
params.dbnsfp_tbi           = getGenomeAttribute('dbnsfp_tbi')
params.spliceai_indel       = getGenomeAttribute('spliceai_indel')
params.spliceai_indel_tbi   = getGenomeAttribute('spliceai_indel_tbi')
params.spliceai_snv         = getGenomeAttribute('spliceai_snv')
params.spliceai_snv_tbi     = getGenomeAttribute('spliceai_snv_tbi')
params.mastermind           = getGenomeAttribute('mastermind')
params.mastermind_tbi       = getGenomeAttribute('mastermind_tbi')
params.eog                  = getGenomeAttribute('eog')
params.eog_tbi              = getGenomeAttribute('eog_tbi')
params.alphamissense        = getGenomeAttribute('alphamissense')
params.alphamissense_tbi    = getGenomeAttribute('alphamissense_tbi')
params.vcfanno_resources    = getGenomeAttribute('vcfanno_resources')
params.vcfanno_config       = getGenomeAttribute('vcfanno_config')
params.hapmap               = getGenomeAttribute('hapmap')           
params.hapmap_tbi           = getGenomeAttribute('hapmap_tbi')
params.omni_1000G            = getGenomeAttribute('omni_1000G')
params.omni_1000G_tbi        = getGenomeAttribute('omni_1000G_tbi')
params.snps_1000G            = getGenomeAttribute('snps_1000G')
params.snps_1000G_tbi        = getGenomeAttribute('snps_1000G_tbi')
params.indels_1000G          = getGenomeAttribute('indels_1000G')
params.indels_1000G_tbi      = getGenomeAttribute('indels_1000G_tbi')
params.known_indels         = getGenomeAttribute('known_indels')
params.known_indels_tbi     = getGenomeAttribute('known_indels_tbi')

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
for(caller in callers) {
    if(!(caller in GlobalVariables.availableCallers)) { error("\"${caller}\" is not a supported callers please use one or more of these instead: ${GlobalVariables.availableCallers}")}
}

if (params.output_suffix && callers.size() > 1) {
    error("Cannot use --output_suffix with more than one caller")
}

if(params.vqsr) {
    if(!params.hapmap || !params.hapmap_tbi) {
        error("Please provide --hapmap and --hapmap_tbi when using --vqsr")
    }
    if(!params.omni_1000G || !params.omni_1000G_tbi) {
        error("Please provide --omni_1000G and --omni_1000G_tbi when using --vqsr")
    }
    if(!params.snps_1000G || !params.snps_1000G_tbi) {
        error("Please provide --snps_1000G and --snps_1000G_tbi when using --vqsr")
    }
    if(!params.indels_1000G || !params.indels_1000G_tbi) {
        error("Please provide --indels_1000G and --indels_1000G_tbi when using --vqsr")
    }
    if(!params.known_indels || !params.known_indels_tbi) {
        error("Please provide --known_indels and --known_indels_tbi when using --vqsr")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

multiqc_logo     = params.multiqc_logo   ?: "$projectDir/assets/CMGG_logo.png"


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
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GERMLINE (
        // Input channels
        samplesheet,

        // File inputs
        params.fasta,
        params.fai,
        params.dict,
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
        params.hapmap,
        params.hapmap_tbi,
        params.omni_1000G,
        params.omni_1000G_tbi,
        params.snps_1000G,
        params.snps_1000G_tbi,
        params.indels_1000G,
        params.indels_1000G_tbi,
        params.known_indels,
        params.known_indels_tbi,

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
        params.updio,
        params.automap,
        params.vep_dbnsfp,
        params.vep_spliceai,
        params.vep_mastermind,
        params.vep_eog,
        params.vep_alphamissense,
        params.vqsr,

        // Value inputs
        params.genome,
        params.species,
        params.vep_cache_version,
        params.vep_chunk_size,
        params.scatter_count,
        params.callers.tokenize(",")
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
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCMGG_GERMLINE (
        PIPELINE_INITIALISATION.out.samplesheet
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
