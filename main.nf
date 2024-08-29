#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/germline
----------------------------------------------------------------------------------------
*/

// Enables the workflow output definition: https://www.nextflow.io/docs/latest/workflow.html#workflow-output-def
nextflow.preview.output = true

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
include { samplesheetToList       } from 'plugin/nf-schema'

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
    vcf_tbi             = GERMLINE.out.vcf_tbi        // channel: [ val(meta), path(vcf), path(tbi) ]
    multiqc_report      = GERMLINE.out.multiqc_report // channel: /path/to/multiqc/report.html
    validation          = GERMLINE.out.validation
    individual_reports  = GERMLINE.out.individual_reports
    family_reports      = GERMLINE.out.family_reports
    individuals_bed     = GERMLINE.out.individuals_bed
    family_bed          = GERMLINE.out.family_bed
    gvcf_tbi            = GERMLINE.out.gvcf_tbi
    updio               = GERMLINE.out.updio
    automap             = GERMLINE.out.automap
    db                  = GERMLINE.out.db
    ped                 = GERMLINE.out.ped

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

    // TODO: remove this once dynamic publish paths have been added to nextflow
    workflow.onComplete = {
        def date = params.skip_date_project ? "" : "${new Date().format("yyyy-MM-dd")}_"
        def final_output = "${params.outdir}/${params.project ? "${date}${params.project}" : "${date}${workflow.runName}"}"
        def ids = samplesheetToList(params.input, "assets/schema_input.json").collect { entry ->
                [ entry[0].id, entry[0].family ]
            }
            .flatten()
            .findAll { id -> id instanceof String && id.length() > 0 }
            .unique()

        // Move around the output directory
        file(params.outdir).eachFileRecurse { file ->
            if (file.isDirectory()) {
                return
            }
            def file_name = file.name
            def file_full_name = file.toString()
            def caller = file_full_name.contains("haplotypecaller") ? "haplotypecaller" :
                file_full_name.contains("vardict") ? "vardict" : ""
            def id = ids.find { id_ss -> file_name.contains(id_ss) } ?: ""
            def extension = file_name.replace("${id}.", "").replace("${caller}.", "")
            if (file_full_name.contains("/temp/vcfs/")) {
                file.moveTo("${final_output}/${id}/${id}.${caller}.${extension}")
            }
            else if (file_full_name.contains("/temp/validation/")) {
                def validation_file = file_name.replace("${caller}.", "")
                file.moveTo("${params.outdir}/${id}/validation/${caller}/${validation_file}")
            }
            else if (file_full_name.contains("/temp/individuals_reports/")) {
                file.moveTo("${params.outdir}/${id}/reports/${file_name}")
            }
            else if (file_full_name.contains("/temp/family_reports/")) {
                file.moveTo("${final_output}/${id}/reports/${file_name}")
            }
            else if (file_full_name.contains("/temp/individuals_beds/")) {
                file.moveTo("${params.outdir}/${id}/${id}.bed")
            }
            else if (file_full_name.contains("/temp/family_beds/")) {
                file.moveTo("${final_output}/${id}/${id}.bed")
            }
            else if (file_full_name.contains("/temp/gvcfs/")) {
                file.moveTo("${params.outdir}/${id}/${id}.${caller}.${extension}")
            }
            else if (file_full_name.contains("/temp/updio/")) {
                def sample = id
                id = file_full_name.split("/temp/updio/")[-1].split("/")[0].replace("updio_${caller}_", "")
                file.moveTo("${final_output}/${id}/updio_${caller}/${sample}/${file_name}")
            }
            else if (file_full_name.contains("/temp/automap/")) {
                def sample = id
                id = file_full_name.split("/temp/automap/")[-1].split("/")[0].replace("automap_${caller}_", "")
                file.moveTo("${final_output}/${id}/automap_${caller}/${sample}/${file_name}")
            }
            else if (file_full_name.contains("/temp/ped/")) {
                file.moveTo("${final_output}/${id}/${id}.${caller}.ped")
            }
            else if (file_full_name.contains("/temp/db/")) {
                file.moveTo("${final_output}/${id}/${id}.${caller}.db")
            }
        }
        file("${params.outdir}/temp").deleteDir()
    }

    publish:
    NFCMGG_GERMLINE.out.vcf_tbi             >> 'temp/vcfs/'
    NFCMGG_GERMLINE.out.validation          >> 'temp/validation/'
    NFCMGG_GERMLINE.out.individual_reports  >> 'temp/individuals_reports/'
    NFCMGG_GERMLINE.out.family_reports      >> 'temp/family_reports/'
    NFCMGG_GERMLINE.out.individuals_bed     >> 'temp/individuals_beds/'
    NFCMGG_GERMLINE.out.family_bed          >> 'temp/family_beds/'
    NFCMGG_GERMLINE.out.gvcf_tbi            >> 'temp/gvcfs/'
    NFCMGG_GERMLINE.out.multiqc_report      >> 'multiqc/'
    NFCMGG_GERMLINE.out.updio               >> 'temp/updio/'
    NFCMGG_GERMLINE.out.automap             >> 'temp/automap/'
    NFCMGG_GERMLINE.out.db                  >> 'temp/db/'
    NFCMGG_GERMLINE.out.ped                 >> 'temp/ped/'
}

output {
    directory "${params.outdir}"
    // TODO: add index once dynamic publish paths have been added to nextflow
    // index {
    //     path 'index.csv'
    // }
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
