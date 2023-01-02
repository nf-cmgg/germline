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
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta                = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fai                  = WorkflowMain.getGenomeAttribute(params, 'fai')
params.dict                 = WorkflowMain.getGenomeAttribute(params, 'dict')
params.strtablefile         = WorkflowMain.getGenomeAttribute(params, 'strtablefile')
params.dbsnp                = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.vep_cache            = WorkflowMain.getGenomeAttribute(params, 'vep_cache')
params.dbnsfp               = WorkflowMain.getGenomeAttribute(params, 'dbnsfp')
params.dbnsfp_tbi           = WorkflowMain.getGenomeAttribute(params, 'dbnsfp_tbi')
params.spliceai_indel       = WorkflowMain.getGenomeAttribute(params, 'spliceai_indel')
params.spliceai_indel_tbi   = WorkflowMain.getGenomeAttribute(params, 'spliceai_indel_tbi')
params.spliceai_snv         = WorkflowMain.getGenomeAttribute(params, 'spliceai_snv')
params.spliceai_snv_tbi     = WorkflowMain.getGenomeAttribute(params, 'spliceai_snv_tbi')
params.mastermind           = WorkflowMain.getGenomeAttribute(params, 'mastermind')
params.mastermind_tbi       = WorkflowMain.getGenomeAttribute(params, 'mastermind_tbi')
params.eog                  = WorkflowMain.getGenomeAttribute(params, 'eog')
params.eog_tbi              = WorkflowMain.getGenomeAttribute(params, 'eog_tbi')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CMGGGERMLINE } from './workflows/cmgg-germline'

//
// WORKFLOW: Run main nf-cmgg-germline analysis pipeline
//
workflow CMGG_CMGGGERMLINE {
    CMGGGERMLINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    CMGG_CMGGGERMLINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
