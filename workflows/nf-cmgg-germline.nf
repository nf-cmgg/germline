/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNfCmggGermline.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input, 
    params.fasta, 
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
    params.vep_merged_cache
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
include { POST_PROCESS             } from '../subworkflows/local/postprocess'
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
include { VCF2DB                                                     } from '../modules/nf-core/modules/vcf2db/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NF_CMGG_GERMLINE {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Importing the pipeline parameters
    //

    fasta              = params.fasta
    genome             = params.genome
    output_mode        = params.output_mode
    species            = params.species
    scatter_count      = params.scatter_count
    use_dragstr_model  = params.use_dragstr_model
    cram_merge         = params.cram_merge

    //
    // Importing the annotation parameters
    //

    vep_cache_version  = params.vep_cache_version
    vep_merged_cache   = params.vep_merged_cache ? params.vep_merged_cache : []

    vep_dbnsfp         = params.vep_dbnsfp
    vep_spliceai       = params.vep_spliceai
    vep_spliceregion   = params.vep_spliceregion
    vep_mastermind     = params.vep_mastermind
    vep_eog            = params.vep_eog

    dbnsfp             = params.dbnsfp
    dbnsfp_tbi         = params.dbnsfp_tbi

    spliceai_snv       = params.spliceai_snv
    spliceai_snv_tbi   = params.spliceai_snv_tbi
    spliceai_indel     = params.spliceai_indel
    spliceai_indel_tbi = params.spliceai_indel_tbi

    mastermind         = params.mastermind
    mastermind_tbi     = params.mastermind_tbi

    eog                = params.eog
    eog_tbi            = params.eog_tbi

    vcfanno            = params.vcfanno
    vcfanno_toml       = params.vcfanno_toml
    vcfanno_resources  = params.vcfanno_resources

    //
    // Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        ch_input
    )

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    inputs = INPUT_CHECK.out.crams
             .multiMap({meta, cram, crai, bed, ped ->
                 new_meta_ped = [:]
                 new_meta_ped.id = meta.family
                 new_meta_ped.family = meta.family

                 new_meta = meta.clone()
                 new_meta.samplename = meta.id

                 beds:                                [new_meta, bed]
                 germline_variant_calling_input_cram: [new_meta, cram, crai]
                 peds:                                [new_meta_ped, ped]
             })

    peds = inputs.peds.distinct()

    //
    // Create the FASTA index from the FASTA file
    //

    if (!params.fasta_fai) {
        FAIDX(
            fasta
        )

        fasta_fai = FAIDX.out.fai
        ch_versions = ch_versions.mix(FAIDX.out.versions)
    } 
    else {
        fasta_fai = params.fasta_fai
    }

    //
    // Create the sequence dictionary from the FASTA file
    //

    if (!params.dict) {
        CREATESEQUENCEDICTIONARY(
            fasta
        )

        dict = CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)
    } 
    else {
        dict = params.dict
    }

    //
    // Create the STR table file from the FASTA file
    //
    if (params.use_dragstr_model) {
        if (!params.strtablefile) {
            COMPOSESTRTABLEFILE(
                fasta,
                fasta_fai,
                dict
            )

            strtablefile = COMPOSESTRTABLEFILE.out.str_table
            ch_versions = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions) 
        } 
        else {
            strtablefile = params.strtablefile
        }
    } 
    else {
        strtablefile = []
    }

    //
    // Perform the variant calling
    //

    GERMLINE_VARIANT_CALLING(
        inputs.germline_variant_calling_input_cram,
        inputs.beds,
        fasta,
        fasta_fai,
        dict,
        strtablefile,
        scatter_count,
        use_dragstr_model,
        cram_merge
    )

    ch_versions = ch_versions.mix(GERMLINE_VARIANT_CALLING.out.versions)

    //
    // Joint-genotyping of the families
    //

    POST_PROCESS(
        GERMLINE_VARIANT_CALLING.out.gvcfs,
        peds,
        fasta,
        fasta_fai,
        dict,
        output_mode
    )

    ch_versions = ch_versions.mix(POST_PROCESS.out.versions)

    //
    // Quality control of the called variants
    //

    VCF_QC(
        POST_PROCESS.out.post_processed_vcfs
    )

    ch_versions = ch_versions.mix(VCF_QC.out.versions)
    ch_reports  = ch_reports.mix(VCF_QC.out.bcftools_stats.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_count.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_qual.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_filter_summary.collect{it[1]}.ifEmpty([]))

    //
    // Annotation of the variants
    //

    if (output_mode == "seqplorer") {
        // Check for the presence of plugins that use extra files
        if (vep_dbnsfp || vep_spliceai || vep_mastermind || vep_eog) {
            vep_extra_files = Channel.empty()
        }
        else {
            vep_extra_files = []
        }

        // Check if all dbnsfp files are given
        if (dbnsfp && dbnsfp_tbi && vep_dbnsfp) {
            vep_extra_files = vep_extra_files.mix(
                Channel.fromPath(params.dbnsfp),
                Channel.fromPath(params.dbnsfp_tbi)
            ).collect()
        }
        else if (dbnsfp || dbnsfp_tbi || vep_dbnsfp) {
            exit 1, "Please specify '--vep_dbsnf true', '--dbnsfp PATH/TO/DBNSFP/FILE' and '--dbnspf_tbi PATH/TO/DBNSFP/INDEX/FILE' to use the dbnsfp VEP plugin."
        }

        // Check if all spliceai files are given
        if (spliceai_snv && spliceai_snv_tbi && spliceai_indel && spliceai_indel_tbi && vep_spliceai) {
            vep_extra_files = vep_extra_files.mix(
                Channel.fromPath(params.spliceai_indel),
                Channel.fromPath(params.spliceai_indel_tbi),
                Channel.fromPath(params.spliceai_snv),
                Channel.fromPath(params.spliceai_snv_tbi)
            ).collect()
        }
        else if (spliceai_snv || spliceai_snv_tbi || spliceai_indel || spliceai_indel_tbi || vep_spliceai) {
            exit 1, "Please specify '--vep_spliceai true', '--spliceai_snv PATH/TO/SPLICEAI/SNV/FILE', '--spliceai_snv_tbi PATH/TO/SPLICEAI/SNV/INDEX/FILE', '--spliceai_indel PATH/TO/SPLICEAI/INDEL/FILE' and '--spliceai_indel_tbi PATH/TO/SPLICEAI/INDEL/INDEX/FILE' to use the SpliceAI VEP plugin."
        }

        // Check if all mastermind files are given
        if (mastermind && mastermind_tbi && vep_mastermind) {
            vep_extra_files = vep_extra_files.mix(
                Channel.fromPath(params.mastermind),
                Channel.fromPath(params.mastermind_tbi)
            ).collect()
        }
        else if (mastermind || mastermind_tbi || vep_mastermind) {
            exit 1, "Please specify '--vep_mastermind true', '--mastermind PATH/TO/MASTERMIND/FILE' and '--mastermind_tbi PATH/TO/MASTERMIND/INDEX/FILE' to use the mastermind VEP plugin."
        }

        // Check if all EOG files are given
        if (eog && eog_tbi && vep_eog) {
            vep_extra_files = vep_extra_files.mix(
                Channel.fromPath(params.eog),
                Channel.fromPath(params.eog_tbi)
            ).collect()
        }
        else if (eog || eog_tbi || vep_eog) {
            exit 1, "Please specify '--vep_eog true', '--eog PATH/TO/EOG/FILE' and '--eog_tbi PATH/TO/EOG/INDEX/FILE' to use the EOG custom VEP plugin."
        }

        // Perform the annotation
        ANNOTATION(
            POST_PROCESS.out.post_processed_vcfs,
            fasta,
            genome,
            species,
            vep_cache_version,
            vep_merged_cache,
            vep_extra_files,
            vcfanno,
            vcfanno_toml,
            vcfanno_resources
        )

        ch_versions = ch_versions.mix(ANNOTATION.out.versions)
        ch_reports  = ch_reports.mix(ANNOTATION.out.reports)
    }  

    //
    // Create Gemini-compatible database files
    //

    if (output_mode == "seqplorer") {
        vcf2db_input = ANNOTATION.out.annotated_vcfs
                       .combine(peds, by: 0)
        
        VCF2DB(
            vcf2db_input
        )
    }

    //
    // Dump the software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_versions_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    //
    // Perform multiQC on all QC data
    //

    ch_multiqc_files = Channel.empty()
    
    ch_multiqc_files = ch_multiqc_files.mix(
                                        ch_versions_yaml,
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
