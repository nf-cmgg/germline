/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

//
// Validate input parameters
//

WorkflowNfCmggGermline.initialise(params, log)

//
// Check input path parameters to see if they exist
//

def checkPathParamList = [ 
    params.fasta, 
    params.fasta_fai,
    params.dict,
    params.strtablefile,
    params.vep_merged_cache,
    params.vcfanno_toml,
    params.vcfanno_resources
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

//
// Check the input samplesheet
//

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }

//
// Check for dependencies between parameters
//

if (params.output_mode == "seqplorer") {
    // Check if a genome is given
    if (!params.genome) { exit 1, "A genome should be supplied for seqplorer mode (use --genome)"}

    // Check if the VEP versions were given
    if (!params.vep_version) { exit 1, "A VEP version should be supplied for seqplorer mode (use --vep_version)"}
    if (!params.vep_cache_version) { exit 1, "A VEP cache version should be supplied for seqplorer mode (use --vep_cache_version)"}

    // Check if a species is entered
    if (!params.species) { exit 1, "A species should be supplied for seqplorer mode (use --species)"}
    
    // Check if all vcfanno files are supplied when vcfanno should be used
    if (params.vcfanno && (!params.vcfanno_toml || !params.vcfanno_resources)) {
        exit 1, "A TOML file and resource directory should be supplied when using vcfanno (use --vcfanno_toml and --vcfanno_resources)"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT THE INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Importing the file pipeline parameters
//

// Input files
fasta              = params.fasta               ? Channel.fromPath(params.fasta).collect()              : Channel.empty()
fasta_fai          = params.fasta_fai           ? Channel.fromPath(params.fasta_fai).collect()          : null
dict               = params.dict                ? Channel.fromPath(params.dict).collect()               : null
strtablefile       = params.strtablefile        ? Channel.fromPath(params.strtablefile).collect()       : null

// Input values
output_mode        = params.output_mode         ?: Channel.empty()
scatter_count      = params.scatter_count       ?: Channel.empty()

// Booleans
always_use_cram    = params.always_use_cram
use_dragstr_model  = params.use_dragstr_model
skip_genotyping    = params.skip_genotyping
use_bcftools_merge = params.use_bcftools_merge

//
// Importing the value pipeline parameters
//

genome             = params.genome              ?: Channel.empty()

//
// Importing the annotation parameters
//

vep_cache_version  = params.vep_cache_version   ?: Channel.empty()
species            = params.species             ?: Channel.empty()

vep_merged_cache   = params.vep_merged_cache    ? Channel.fromPath(params.vep_merged_cache).collect()   : []

vcfanno            = params.vcfanno             ?: Channel.empty()

vcfanno_toml       = params.vcfanno_toml        ? Channel.fromPath(params.vcfanno_toml).collect()       : Channel.empty()
vcfanno_res_inp    = params.vcfanno_resources   ? Channel.fromPath(params.vcfanno_resources).collect()  : Channel.empty()

//
// Check for the presence of EnsemblVEP plugins that use extra files
//

vep_extra_files = []

// Check if all dbnsfp files are given
if (params.dbnsfp && params.dbnsfp_tbi && params.vep_dbnsfp) {
    vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
    vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
}
else if (params.dbnsfp || params.dbnsfp_tbi || params.vep_dbnsfp) {
    exit 1, "Please specify '--vep_dbsnf true', '--dbnsfp PATH/TO/DBNSFP/FILE' and '--dbnspf_tbi PATH/TO/DBNSFP/INDEX/FILE' to use the dbnsfp VEP plugin."
}

// Check if all spliceai files are given
if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi && params.vep_spliceai) {
    vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
}
else if (params.spliceai_snv || params.spliceai_snv_tbi || params.spliceai_indel || params.spliceai_indel_tbi || params.vep_spliceai) {
    exit 1, "Please specify '--vep_spliceai true', '--spliceai_snv PATH/TO/SPLICEAI/SNV/FILE', '--spliceai_snv_tbi PATH/TO/SPLICEAI/SNV/INDEX/FILE', '--spliceai_indel PATH/TO/SPLICEAI/INDEL/FILE' and '--spliceai_indel_tbi PATH/TO/SPLICEAI/INDEL/INDEX/FILE' to use the SpliceAI VEP plugin."
}

// Check if all mastermind files are given
if (params.mastermind && params.mastermind_tbi && params.vep_mastermind) {
    vep_extra_files.add(file(params.mastermind, checkIfExists: true))
    vep_extra_files.add(file(params.mastermind_tbi, checkIfExists: true))
}
else if (params.mastermind || params.mastermind_tbi || params.vep_mastermind) {
    exit 1, "Please specify '--vep_mastermind true', '--mastermind PATH/TO/MASTERMIND/FILE' and '--mastermind_tbi PATH/TO/MASTERMIND/INDEX/FILE' to use the mastermind VEP plugin."
}

// Check if all EOG files are given
if (params.eog && params.eog_tbi && params.vep_eog) {
    vep_extra_files.add(file(params.eog, checkIfExists: true))
    vep_extra_files.add(file(params.eog_tbi, checkIfExists: true))
}
else if (params.eog || params.eog_tbi || params.vep_eog) {
    exit 1, "Please specify '--vep_eog true', '--eog PATH/TO/EOG/FILE' and '--eog_tbi PATH/TO/EOG/INDEX/FILE' to use the EOG custom VEP plugin."
}

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
include { UNTAR                                                      } from '../modules/nf-core/modules/untar/main'
include { VCF2DB                                                     } from '../modules/nf-core/modules/vcf2db/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC                                                    } from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// The main workflow
workflow NF_CMGG_GERMLINE {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Create the optional input files if they are not supplied
    //

    if (!fasta_fai) {
        fasta_fai   = FAIDX(fasta).fai
        ch_versions = ch_versions.mix(FAIDX.out.versions)
    }

    if (!dict) {
        dict        = CREATESEQUENCEDICTIONARY(fasta).dict
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (use_dragstr_model && !strtablefile) {
        strtablefile = COMPOSESTRTABLEFILE(fasta,fasta_fai,dict).str_table
        ch_versions  = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions)
    }

    if (output_mode == "seqplorer" && vcfanno && params.vcfanno_resources.endsWith(".tar.gz")) {
        vcfanno_resources = UNTAR( vcfanno_res_inp.map({dir -> [ [], dir ]}) ).untar.map({meta, dir -> dir})
        ch_versions       = ch_versions.mix(UNTAR.out.versions)
    } else {
        vcfanno_resources = vcfanno_res_inp
    }

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
        always_use_cram
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
        output_mode,
        skip_genotyping,
        use_bcftools_merge
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
