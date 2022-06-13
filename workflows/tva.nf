/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTva.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input, 
    params.fasta, 
    params.dbnsfp, 
    params.dbsnfp_tbi,
    params.spliceai_indel,
    params.spliceai_indel_tbi,
    params.spliceai_snv,
    params.spliceai_snv_tbi,
    params.mastermind,
    params.mastermind_tbi,
    params.eog,
    params.eog_tbi
]
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
include { VCF_QC                   } from '../subworkflows/local/vcf_qc'
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
include { MULTIQC                                                    } from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TVA {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

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
    // Quality control of the called variants
    //

    VCF_QC(
        GENOTYPE.out.genotyped_vcfs
    )

    ch_versions = ch_versions.mix(VCF_QC.out.versions)
    ch_reports  = ch_reports.mix(VCF_QC.out.bcftools_stats.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_count.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_qual.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_filter_summary.collect{it[1]}.ifEmpty([]))

    //
    // Annotation of the variants
    //

    vep_extra_files = Channel.empty()

    if (params.dbnsfp && params.dbnsfp_tbi) {
        vep_extra_files = vep_extra_files.mix(
            Channel.fromPath(params.dbnsfp),
            Channel.fromPath(params.dbnsfp_tbi)
        ).collect()
    }

    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
        vep_extra_files = vep_extra_files.mix(
            Channel.fromPath(params.spliceai_indel),
            Channel.fromPath(params.spliceai_indel_tbi),
            Channel.fromPath(params.spliceai_snv),
            Channel.fromPath(params.spliceai_snv_tbi)
        ).collect()
    }

    if (params.mastermind && params.mastermind_tbi) {
        vep_extra_files = vep_extra_files.mix(
            Channel.fromPath(params.mastermind),
            Channel.fromPath(params.mastermind_tbi)
        ).collect()
    }

    if (params.eog && params.eog_tbi) {
        vep_extra_files = vep_extra_files.mix(
            Channel.fromPath(params.eog),
            Channel.fromPath(params.eog_tbi)
        ).collect()
    }

    ANNOTATION(
        GENOTYPE.out.genotyped_vcfs,
        fasta,
        vep_extra_files
    )

    ch_versions = ch_versions.mix(ANNOTATION.out.versions)
    ch_reports  = ch_reports.mix(ANNOTATION.out.reports)  

    //
    // Dump the software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    //
    // Perform multiQC on all QC data
    //

    ch_multiqc_files = Channel.empty()
    
    ch_multiqc_files = ch_multiqc_files.mix(
                                        ch_version_yaml,
                                        ch_reports.collect(),
                                        ch_multiqc_custom_config
                                        )
    MULTIQC(
        ch_multiqc_files.collect()
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
