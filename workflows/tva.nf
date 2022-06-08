/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTva.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK              } from '../subworkflows/local/input_check'
include { GERMLINE_VARIANT_CALLING } from '../subworkflows/local/germline_variant_calling'
include { GENOTYPE                 } from '../subworkflows/local/genotype'
include { ANNOTATION               } from '../subworkflows/local/annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FAIDX as FAIDX                                    } from '../modules/nf-core/modules/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as CREATESEQUENCEDICTIONARY } from '../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { GATK4_COMPOSESTRTABLEFILE as COMPOSESTRTABLEFILE           } from '../modules/nf-core/modules/gatk4/composestrtablefile/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TVA {

    ch_versions = Channel.empty()

    //
    // Importing the parameters
    //

    fasta = params.fasta

    //
    // Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Create the FASTA index from the FASTA file
    //

    if (!params.fasta_fai){
        FAIDX(
            fasta
        )
        fasta_fai = FAIDX.out.fai
        ch_versions = ch_versions.mix(FAIDX.out.versions)

    } else {
        fasta_fai = params.fasta_fai
    }

    //
    // Create the sequence dictionary from the FASTA file
    //

    if (!params.dict){
        CREATESEQUENCEDICTIONARY(
            fasta
        )
        dict = CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)
    } else {
        dict = params.dict
    }

    //
    // Create the STR table file from the FASTA file
    //
    if (params.use_dragstr_model){
        if (!params.strtablefile){
            COMPOSESTRTABLEFILE(
                fasta,
                fasta_fai,
                dict
            )
            strtablefile = COMPOSESTRTABLEFILE.out.str_table
            ch_versions = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions)
        } else {
            strtablefile = params.strtablefile
        }
    } else {
        strtablefile = []
    }

    //
    // Perform the variant calling
    //

    beds = INPUT_CHECK.out.crams.map(
    {meta, cram, crai, bed ->
        [meta, bed]
    })

    germline_variant_calling_input_cram = INPUT_CHECK.out.crams.map(
    {meta, cram, crai, bed ->
        [meta, cram, crai]
    })

    GERMLINE_VARIANT_CALLING(
        germline_variant_calling_input_cram,
        beds,
        fasta,
        fasta_fai,
        dict,
        strtablefile
    )

    ch_versions = ch_versions.mix(GERMLINE_VARIANT_CALLING.out.versions)

    //
    // Joint-genotyping of the families
    //

    GENOTYPE(
        GERMLINE_VARIANT_CALLING.out.vcfs,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(GENOTYPE.out.versions)

    //
    // Annotation of the variants
    //

    ANNOTATION(
        GENOTYPE.out.genotyped_gvcfs
    )

    ch_versions = ch_versions.mix(ANNOTATION.out.versions)    

    //
    // Dump the software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
