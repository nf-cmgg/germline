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

vcfanno            = params.vcfanno

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

ch_multiqc_config   = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
multiqc_logo        = params.multiqc_logo   ? file(params.multiqc_logo, checkIfExists: true)   : file("$projectDir/assets/CMGG_logo.png", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

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

include { SAMTOOLS_FAIDX as FAIDX                                    } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_COMPOSESTRTABLEFILE as COMPOSESTRTABLEFILE           } from '../modules/nf-core/gatk4/composestrtablefile/main'
include { INDEX_TO_BED                                               } from '../modules/local/indextobed'
include { UNTAR                                                      } from '../modules/nf-core/untar/main'
include { VCF2DB                                                     } from '../modules/nf-core/vcf2db/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                                                    } from '../modules/nf-core/multiqc/main'

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
        FAIDX(
            fasta.map({ fasta -> [ [id:"fasta_fai"], fasta ]})
        )
        ch_versions = ch_versions.mix(FAIDX.out.versions)

        FAIDX.out.fai
            .map({ meta, fasta_fai -> [ fasta_fai ]})
            .collect()
            .set { fasta_fai }
    }

    if (!dict) {
        CREATESEQUENCEDICTIONARY(
            fasta
        )
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)

        CREATESEQUENCEDICTIONARY.out.dict
            .collect()
            .set { dict }
    }

    if (use_dragstr_model && !strtablefile) {
        COMPOSESTRTABLEFILE(
            fasta,
            fasta_fai,
            dict
        )
        ch_versions  = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions)

        COMPOSESTRTABLEFILE.out.str_table
            .collect()
            .set { strtablefile }
    }

    if (output_mode == "seqplorer" && vcfanno) {
        if (params.vcfanno_resources.endsWith(".tar.gz")) {
            UNTAR(
                vcfanno_res_inp.map({dir -> [ [], dir ]})
            )
            ch_versions = ch_versions.mix(UNTAR.out.versions)

            UNTAR.out.untar
                .map(
                    { meta, dir ->
                        dir
                    }
                )
                .collect()
                .set { vcfanno_resources }
        } else {
            Channel.value(vcfanno_res_inp)
                .set { vcfanno_resources }
        }
    } else {
        vcfanno_resources = []
    }

    //
    // Read in samplesheet, validate and stage input files
    //

    parse_input(ch_input)
        .multiMap(
            { meta, cram, crai, bed, ped ->
                ped_family_id = meta.family ?: get_family_id_from_ped(ped)

                new_meta_ped = [:]
                new_meta = meta.clone()

                new_meta_ped.id     = meta.family ?: ped_family_id
                new_meta_ped.family = meta.family ?: ped_family_id
                new_meta.family     = meta.family ?: ped_family_id

                beds:                                [new_meta, bed]
                germline_variant_calling_input_cram: [new_meta, cram, crai]
                peds:                                [new_meta_ped, ped]
            }
        )
        .set { ch_parsed_inputs }

    // TODO specify the family size (warning: a sample can be supplied multiple times. This should not be counted towards the family total)

    ch_parsed_inputs.beds
        .branch(
            { meta, bed ->
                valid: bed != []
                invalid: bed == []
            }
        )
        .set { beds }

    //
    // Create a BED file from the FASTA index to paralellize genome calls
    //

    // TODO Only perform this once

    INDEX_TO_BED(
        beds.invalid.combine(fasta_fai).map({ meta, bed, fasta_fai -> [ meta, fasta_fai ]})
    )

    ch_parsed_inputs.peds
        .distinct()
        .groupTuple()
        .map(
            { meta, peds ->
                output = peds.size() == 1 ? peds[0] : peds[0] == [] ? peds[1] : peds[0]
                [ meta, output ]
            }
        )
        .set { peds }

    //
    // Perform the variant calling
    //

    GERMLINE_VARIANT_CALLING(
        ch_parsed_inputs.germline_variant_calling_input_cram,
        beds.valid.mix(INDEX_TO_BED.out.bed),
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
        ANNOTATION.out.annotated_vcfs
            .combine(peds, by: 0)
            .set { vcf2db_input }

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

    ch_multiqc_files = ch_multiqc_files.mix(ch_versions_yaml, ch_reports.collect())

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        [],
        multiqc_logo
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
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parse_input(input_csv) {

    // The samplesheet schema (change this to adjust the input check)
    def samplesheet_schema = [
        'columns': [
            'sample': [
                'content': 'meta',
                'meta_name': 'id,sample',
                'pattern': '',
            ],
            'family': [
                'content': 'meta',
                'meta_name': 'family',
                'pattern': ''
            ],
            'cram': [
                'content': 'file',
                'pattern': '^.*\\.cram$',
            ],
            'crai': [
                'content': 'file',
                'pattern': '^.*\\.crai$',
            ],
            'bed': [
                'content': 'file',
                'pattern': '^.*\\.bed$',
            ],
            'ped': [
                'content': 'file',
                'pattern': '^.*\\.ped$',
            ]
        ],
        'required': ['sample','cram'],
    ]

    // Don't change these variables
    def row_count = 1
    def all_columns = samplesheet_schema.columns.keySet().collect()
    def mandatory_columns = samplesheet_schema.required

    // Header checks
    Channel.value(input_csv).splitCsv(strip:true).first().map({ row ->

        if(row != all_columns) {
            def commons = all_columns.intersect(row)
            def diffs = all_columns.plus(row)
            diffs.removeAll(commons)

            if(diffs.size() > 0){
                def missing_columns = []
                def wrong_columns = []
                for(diff : diffs){
                    diff in all_columns ? missing_columns.add(diff) : wrong_columns.add(diff)
                }
                if(missing_columns.size() > 0){
                    exit 1, "[Samplesheet Error] The column(s) $missing_columns is/are not present. The header should look like: $all_columns"
                }
                else {
                    exit 1, "[Samplesheet Error] The column(s) $wrong_columns should not be in the header. The header should look like: $all_columns"
                }
            }
            else {
                exit 1, "[Samplesheet Error] The columns $row are not in the right order. The header should look like: $all_columns"
            }

        }
    })

    // Field checks + returning the channels
    Channel.value(input_csv).splitCsv(header:true, strip:true).map({ row ->

        row_count++

        // Check the mandatory columns
        def missing_mandatory_columns = []
        for(column : mandatory_columns) {
            row[column] ?: missing_mandatory_columns.add(column)
        }
        if(missing_mandatory_columns.size > 0){
            exit 1, "[Samplesheet Error] The mandatory column(s) $missing_mandatory_columns is/are empty on line $row_count"
        }

        def output = []
        def meta = [:]
        for(col : samplesheet_schema.columns) {
            key = col.key
            content = row[key]

            if(!(content ==~ col.value['pattern']) && col.value['pattern'] != '' && content != '') {
                exit 1, "[Samplesheet Error] The content of column '$key' on line $row_count does not match the pattern '${col.value['pattern']}'"
            }

            if(col.value['content'] == 'file'){
                output.add(content ? file(content, checkIfExists:true) : [])
            }
            else if(col.value['content'] == 'meta' && content != ''){
                for(meta_name : col.value['meta_name'].split(",")){
                    meta[meta_name] = content.replace(' ', '_')
                }
            }
        }

        output.add(0, meta)
        return output
    })

}

def get_family_id_from_ped(ped_file){

    // Check if there is a file
    if (ped_file.isEmpty()){
        exit 1, "[PED file error] Please specify either the family or enter a PED file into the samplesheet for each sample"
    }

    // Read the PED file
    def ped = file(ped_file, checkIfExists: true).text

    // Perform a validity check on the PED file since vcf2db is picky and not capable of giving good error messages
    comment_count = 0
    line_count = 0

    for( line : ped.readLines()) {
        line_count++
        if (line_count == 1 && line ==~ /^#.*$/) {
            continue
        }
        else if (line_count > 1 && line ==~ /^#.*$/) {
            exit 1, "[PED file error] A commented line was found on line ${line_count} in ${ped_file}, the only commented line allowed is an optional header on line 1."
        }
        else if (line_count == 1 && line ==~ /^#.* $/) {
            exit 1, "[PED file error] The header in ${ped_file} contains a trailing space, please remove this."
        }
        else if (line ==~ /^.+#.*$/) {
            exit 1, "[PED file error] A '#' has been found as a non-starting character on line ${line_count} in ${ped_file}, this is an illegal character and should be removed."
        }
        else if (line ==~ /^[^#].* .*$/) {
            exit 1, "[PED file error] A space has been found on line ${line_count} in ${ped_file}, please only use tabs to seperate the values (and change spaces in names to '_')."
        }
        else if ((line ==~ /^(\w+\t)+\w+$/) == false) {
            exit 1, "[PED file error] An illegal character has been found on line ${line_count} in ${ped_file}, only a-z; A-Z; 0-9 and '_' are allowed as column values."
        }
        else if ((line ==~ /^(\w+\t){5}\w+$/) == false) {
            exit 1, "[PED file error] ${ped_file} should contain exactly 6 tab-delimited columns (family_id    individual_id    paternal_id    maternal_id    sex    phenotype). This is not the case on line ${line_count}."
        }
    }
    if (ped =~ /\n$/) {
        exit 1, "[PED file error] An empty new line has been detected at the end of ${ped_file}, please remove this line."
    }

    // get family_id
    return (ped =~ /\n([^#]\w+)/)[0][1]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
