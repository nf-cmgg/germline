#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CenterForMedicalGeneticsGhent/nf-cmgg-germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-germline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CMGGGERMLINE            } from './workflows/cmgg-germline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_cmgg_germline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_cmgg_germline_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_cmgg_germline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CMGG_CMGGGERMLINE {

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
    CMGGGERMLINE (
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
        CMGGGERMLINE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    CMGG_CMGGGERMLINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
