/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

//
// Validate input parameters
//

WorkflowCmggGermline.initialise(params, log)

//
// Check input path parameters to see if they exist
//

def checkPathParamList = [
    params.fasta,
    params.fai,
    params.dict,
    params.strtablefile,
    params.vep_cache,
    params.vcfanno_config,
    params.vcfanno_lua,
    params.dbsnp,
    params.dbsnp_tbi,
    params.somalier_sites,
    params.sdf
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

//
// Check the input samplesheet
//

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }

//
// Check for dependencies between parameters
//

if(params.dbsnp_tbi && !params.dbsnp){
    exit 1, "Please specify the dbsnp VCF with --dbsnp VCF"
}

if (params.annotate) {
    // Check if a genome is given
    if (!params.genome) { exit 1, "A genome should be supplied for seqplorer mode (use --genome)"}

    // Check if the VEP versions were given
    if (!params.vep_version) { exit 1, "A VEP version should be supplied for seqplorer mode (use --vep_version)"}
    if (!params.vep_cache_version) { exit 1, "A VEP cache version should be supplied for seqplorer mode (use --vep_cache_version)"}

    // Check if a species is entered
    if (!params.species) { exit 1, "A species should be supplied for seqplorer mode (use --species)"}

    // Check if all vcfanno files are supplied when vcfanno should be used
    if (params.vcfanno && (!params.vcfanno_config || !params.vcfanno_resources)) {
        exit 1, "A TOML file and resource files should be supplied when using vcfanno (use --vcfanno_config and --vcfanno_resources)"
    }
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

include { PREPROCESSING                 } from '../subworkflows/local/preprocessing'
include { GERMLINE_VARIANT_CALLING      } from '../subworkflows/local/germline_variant_calling'
include { JOINT_GENOTYPING              } from '../subworkflows/local/joint_genotyping'
include { ANNOTATION                    } from '../subworkflows/local/annotation'
include { ADD_PED_HEADER                } from '../subworkflows/local/add_ped_header'
include { VCF_VALIDATE_SMALL_VARIANTS   } from '../subworkflows/local/vcf_validate_small_variants/main'

include { VCF_EXTRACT_RELATE_SOMALIER   } from '../subworkflows/nf-core/vcf_extract_relate_somalier/main'

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
include { RTGTOOLS_FORMAT                                            } from '../modules/nf-core/rtgtools/format/main'
include { UNTAR                                                      } from '../modules/nf-core/untar/main'
include { TABIX_TABIX as TABIX_DBSNP                                 } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TRUTH                                 } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER as FILTER_SNPS                             } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_INDELS                           } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_FAMILY                    } from '../modules/nf-core/bcftools/stats/main'
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
workflow CMGGGERMLINE {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Importing the input files
    //
    fasta              = Channel.fromPath(params.fasta).collect()
    fasta_fai          = params.fai                 ? Channel.fromPath(params.fai).collect()                                        : null
    dict               = params.dict                ? Channel.fromPath(params.dict).collect()                                       : null
    strtablefile       = params.strtablefile        ? Channel.fromPath(params.strtablefile).collect()                               : null
    sdf                = params.sdf                 ? Channel.fromPath(params.sdf).map {sdf -> [[id:'reference'], sdf]}.collect()   : null

    default_roi        = params.roi                 ? Channel.fromPath(params.roi).collect()                : []

    dbsnp              = params.dbsnp               ? Channel.fromPath(params.dbsnp).collect()              : []
    dbsnp_tbi          = params.dbsnp_tbi           ? Channel.fromPath(params.dbsnp_tbi).collect()          : []

    somalier_sites     = params.somalier_sites      ? Channel.fromPath(params.somalier_sites).collect()     : []

    vep_cache          = params.vep_cache           ? Channel.fromPath(params.vep_cache).collect()          : []

    vcfanno_config     = params.vcfanno_config      ? Channel.fromPath(params.vcfanno_config).collect()     : []
    vcfanno_lua        = params.vcfanno_lua         ? Channel.fromPath(params.vcfanno_lua).collect()        : []
    vcfanno_resources  = params.vcfanno_resources   ? Channel.of(params.vcfanno_resources.split(",")).map({ file(it, checkIfExists:true) }).collect()   : []

    //
    // Check for the presence of EnsemblVEP plugins that use extra files
    //

    if(params.annotate){
        vep_extra_files = []

        // Check if all dbnsfp files are given
        if (params.dbnsfp && params.dbnsfp_tbi && params.vep_dbnsfp) {
            vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
            vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
        }
        else if (params.vep_dbnsfp) {
            exit 1, "Please specify '--vep_dbsnfp true', '--dbnsfp PATH/TO/DBNSFP/FILE' and '--dbnspf_tbi PATH/TO/DBNSFP/INDEX/FILE' to use the dbnsfp VEP plugin."
        }

        // Check if all spliceai files are given
        if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi && params.vep_spliceai) {
            vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
            vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
            vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
            vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        }
        else if (params.vep_spliceai) {
            exit 1, "Please specify '--vep_spliceai true', '--spliceai_snv PATH/TO/SPLICEAI/SNV/FILE', '--spliceai_snv_tbi PATH/TO/SPLICEAI/SNV/INDEX/FILE', '--spliceai_indel PATH/TO/SPLICEAI/INDEL/FILE' and '--spliceai_indel_tbi PATH/TO/SPLICEAI/INDEL/INDEX/FILE' to use the SpliceAI VEP plugin."
        }

        // Check if all mastermind files are given
        if (params.mastermind && params.mastermind_tbi && params.vep_mastermind) {
            vep_extra_files.add(file(params.mastermind, checkIfExists: true))
            vep_extra_files.add(file(params.mastermind_tbi, checkIfExists: true))
        }
        else if (params.vep_mastermind) {
            exit 1, "Please specify '--vep_mastermind true', '--mastermind PATH/TO/MASTERMIND/FILE' and '--mastermind_tbi PATH/TO/MASTERMIND/INDEX/FILE' to use the mastermind VEP plugin."
        }

        // Check if all maxentscan files are given
        if (params.maxentscan && params.vep_maxentscan) {
            vep_extra_files.add(file(params.maxentscan, checkIfExists: true))
        }
        else if (params.vep_maxentscan) {
            exit 1, "Please specify '--vep_maxentscan true', '--maxentscan PATH/TO/MAXENTSCAN/' to use the MaxEntScan VEP plugin."
        }

        // Check if all EOG files are given
        if (params.eog && params.eog_tbi && params.vep_eog) {
            vep_extra_files.add(file(params.eog, checkIfExists: true))
            vep_extra_files.add(file(params.eog_tbi, checkIfExists: true))
        }
        else if (params.vep_eog) {
            exit 1, "Please specify '--vep_eog true', '--eog PATH/TO/EOG/FILE' and '--eog_tbi PATH/TO/EOG/INDEX/FILE' to use the EOG custom VEP plugin."
        }
    }

    //
    // Create the optional input files if they are not supplied
    //

    if (dbsnp != [] && !dbsnp_tbi == []) {
        TABIX_DBSNP(
            dbsnp.map { [[id:'dbsnp'], it] }
        )

        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
        TABIX_DBSNP.out.tbi
            .map(
                { meta, tbi ->
                    [ tbi ]
                }
            )
            .collect()
            .set { dbsnp_tbi }
    }

    if (!fasta_fai) {
        FAIDX(
            fasta.map({ fasta -> [ [id:"fasta_fai"], fasta ]})
        )
        ch_versions = ch_versions.mix(FAIDX.out.versions)

        FAIDX.out.fai
            .map({ meta, fasta_fai -> [ fasta_fai ]})
            .collect()
            .dump(tag:'fasta_fai', pretty:true)
            .set { fasta_fai }
    }

    if (!dict) {
        CREATESEQUENCEDICTIONARY(
            fasta
        )
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)

        CREATESEQUENCEDICTIONARY.out.dict
            .collect()
            .dump(tag:'dict', pretty:true)
            .set { dict }
    }

    if (params.use_dragstr_model && !strtablefile) {
        COMPOSESTRTABLEFILE(
            fasta,
            fasta_fai,
            dict
        )
        ch_versions  = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions)

        COMPOSESTRTABLEFILE.out.str_table
            .collect()
            .dump(tag:'strtablefile', pretty:true)
            .set { strtablefile }
    }

    if (params.validate && !sdf) {
        RTGTOOLS_FORMAT(
            fasta.map { fasta -> [[id:'reference'], fasta, [], []]}
        )
        ch_versions  = ch_versions.mix(RTGTOOLS_FORMAT.out.versions)

        RTGTOOLS_FORMAT.out.sdf
            .collect()
            .dump(tag:'sdf', pretty:true)
            .set { sdf }
    }
    else if (params.validate) {
        sdf.branch { meta, sdf ->
            zip = sdf.name.endsWith(".tar.gz")
            tarzipped: zip
            untarred: !zip
        }
        .set { sdf_branch }

        UNTAR(
            sdf_branch.tarzipped
        )
        ch_versions = ch_versions.mix(UNTAR.out.versions)

        UNTAR.out.untar
            .mix(sdf_branch.untarred)
            .dump(tag:'sdf', pretty:true)
            .set { sdf }
    }

    //
    // Read in samplesheet, validate and convert to a channel
    //

    SamplesheetConversion.convert(ch_input, file("${projectDir}/assets/schema_input.json", checkIfExists:true))
        .map { meta, cram, crai, callable, roi, ped, truth_vcf, truth_tbi ->
            // Infer the family ID from the PED file if no family ID was given, if no PED is given use the sample ID as family ID
            // Infer the type of data (WES or WGS). When a ROI is given through the params or samplesheet, the sample is marked as WES, otherwise it is WGS
            new_meta = meta + [
                family: meta.family ?: ped ? get_family_id_from_ped(ped) : meta.sample, 
            ]
            [ new_meta, cram, crai, callable, roi, ped, truth_vcf, truth_tbi ]
        }
        .tap { ch_raw_inputs }
        .map { it[0] }
        .distinct()
        .map { meta ->
            [ meta.family, 1 ]
        }
        .groupTuple()
        .map { family, count ->
            counts = family ? count.sum() : 1
            [ family, counts ]
        }
        .combine(
            ch_raw_inputs
                .map { meta, cram, crai, callable, roi, ped, truth_vcf, truth_tbi ->
                    [ meta.family, meta, cram, crai, callable, roi, ped, truth_vcf, truth_tbi ]
                }
        , by:0)
        .multiMap(
            { family, family_count, meta, cram, crai, callable, roi, ped, truth_vcf, truth_tbi ->
                new_meta_ped = [
                    id:             meta.family,
                    family:         meta.family,
                    family_count:   family_count
                ]

                new_meta = meta + [family_count:family_count]

                truth_variants: [new_meta, truth_vcf, truth_tbi]
                cram:           [new_meta, cram, crai]
                peds:           [new_meta_ped, ped]
                roi:            [new_meta, roi]
                callable:       [new_meta, callable]
            }
        )
        .set { ch_parsed_inputs }

    ch_parsed_inputs.roi.dump(tag:'input_roi', pretty:true)
    ch_parsed_inputs.callable.dump(tag:'input_callable', pretty:true)
    ch_parsed_inputs.truth_variants.dump(tag:'truth_variants', pretty:true)
    ch_parsed_inputs.cram.dump(tag:'input_crams', pretty:true)
    ch_parsed_inputs.peds.dump(tag:'input_peds', pretty:true)

    //
    // Run sample preprocessing
    //

    PREPROCESSING(
        ch_parsed_inputs.cram,
        ch_parsed_inputs.roi,
        ch_parsed_inputs.callable,
        fasta,
        fasta_fai,
        default_roi
    )

    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    ch_parsed_inputs.peds
        .distinct()
        .groupTuple()
        .map(
            { meta, peds ->
                output = peds.size() == 1 ? peds[0] : peds[0] == [] ? peds[1] : peds[0]
                [ meta, output ]
            }
        )
        .dump(tag:'peds', pretty:true)
        .set { peds }

    //
    // Perform the variant calling
    //

    GERMLINE_VARIANT_CALLING(
        PREPROCESSING.out.ready_crams,
        PREPROCESSING.out.ready_beds,
        fasta,
        fasta_fai,
        dict,
        strtablefile,
        dbsnp,
        dbsnp_tbi
    )

    ch_versions = ch_versions.mix(GERMLINE_VARIANT_CALLING.out.versions)
    ch_reports  = ch_reports.mix(GERMLINE_VARIANT_CALLING.out.reports)

    GERMLINE_VARIANT_CALLING.out.gvcfs
        .dump(tag:'variantcalling_output', pretty:true)
        .set { variantcalling_output }

    //
    // Validate the found variants
    //

    if (params.validate){

        variantcalling_output
            .join(PREPROCESSING.out.ready_beds, failOnDuplicate: true, failOnMismatch: true)
            .join(ch_parsed_inputs.truth_variants, failOnDuplicate: true, failOnMismatch: true)
            .filter { meta, vcf, tbi, bed, truth_vcf, truth_tbi ->
                truth_vcf != []
            }
            .multiMap { meta, vcf, tbi, bed, truth_vcf, truth_tbi ->
                vcf: [ meta, vcf, tbi, truth_vcf, truth_tbi ]
                bed: [ meta, [], bed ]
            }
            .set { validation_input }

        validation_input.vcf
            .branch { meta, vcf, tbi, truth_vcf, truth_tbi ->
                tbi: truth_tbi != []
                no_tbi: truth_tbi == []
            }
            .set { validation_branch }

        TABIX_TRUTH(
            validation_branch.no_tbi.map { meta, vcf, tbi, truth_vcf, truth_tbi -> [ meta, truth_vcf ]}
        )

        ch_versions = ch_versions.mix(TABIX_TRUTH.out.versions)

        validation_branch.no_tbi
            .join(TABIX_TRUTH.out.tbi, failOnDuplicate: true, failOnMismatch: true) 
            .map { meta, vcf, tbi, truth_vcf, empty_tbi, truth_tbi ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi ]
            }
            .mix(validation_branch.tbi)
            .set { validation_vcfs }

        // TODO: add support for truth regions when happy works
        VCF_VALIDATE_SMALL_VARIANTS(
            validation_vcfs,
            validation_input.bed,
            fasta.map { [[], it] },
            fasta_fai.map { [[], it] },
            sdf,
            [[],[]],
            [[],[]],
            [[],[]],
            "vcfeval" //Only VCFeval for now, awaiting the conda fix for happy (https://github.com/bioconda/bioconda-recipes/pull/39267)
        )

        ch_versions = ch_versions.mix(VCF_VALIDATE_SMALL_VARIANTS.out.versions)
    }

    //
    // Joint-genotyping of the families
    //

    JOINT_GENOTYPING(
        variantcalling_output,
        PREPROCESSING.out.ready_crams,
        PREPROCESSING.out.ready_beds,
        peds,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(JOINT_GENOTYPING.out.versions)

    JOINT_GENOTYPING.out.genotyped_vcfs
        .dump(tag:'joint_genotyping_output', pretty:true)
        .tap { vcf_qc_input }
        .set { joint_genotyping_output }

    //
    // Filter the variants
    //

    if (params.filter) {
        FILTER_SNPS(
            joint_genotyping_output
        )

        FILTER_INDELS(
            FILTER_SNPS.out.vcf
        )

        ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)
        ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

        FILTER_INDELS.out.vcf.set { filter_output }

    }
    else {
        joint_genotyping_output.set { filter_output }
    }

    filter_output.dump(tag:'filter_output', pretty: true)

    //
    // Somalier
    //

    VCF_EXTRACT_RELATE_SOMALIER(
        filter_output
            .filter { meta, vcf ->
                meta.family_count > 1
            }
            .map { it + [[], 1] },
        fasta,
        fasta_fai,
        somalier_sites,
        peds
            .filter { meta, ped ->
                meta.family_count > 1
            },
        [],
        []
    )

    ch_versions = ch_versions.mix(VCF_EXTRACT_RELATE_SOMALIER.out.versions)

    //
    // Add PED headers to the VCFs
    //

    if(params.add_ped){
        ADD_PED_HEADER(
            filter_output,
            VCF_EXTRACT_RELATE_SOMALIER.out.samples_tsv
        )

        ch_versions = ch_versions.mix(ADD_PED_HEADER.out.versions)

        ADD_PED_HEADER.out.ped_vcfs
            .dump(tag:'ped_vcfs', pretty:true)
            .set { ped_vcfs }
    } else {
        filter_output.set { ped_vcfs }
    }



    //
    // Annotation of the variants and creation of Gemini-compatible database files
    //
    if (params.annotate) {

        ANNOTATION(
            ped_vcfs,
            fasta,
            fasta_fai,
            vep_cache,
            vep_extra_files,
            vcfanno_config,
            vcfanno_lua,
            vcfanno_resources
        )
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)
        ch_reports  = ch_reports.mix(ANNOTATION.out.reports)

        ANNOTATION.out.annotated_vcfs.set { annotation_output }
    } else {
        ped_vcfs.set { annotation_output }
    }

    annotation_output
        .dump(tag:'annotation_output', pretty:true)
        .set { annotation_output }

    //
    // Perform QC on the final VCFs
    //

    BCFTOOLS_STATS_FAMILY(
        annotation_output.map{ it + [[]] },
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS_FAMILY.out.versions)
    ch_reports  = ch_reports.mix(BCFTOOLS_STATS_FAMILY.out.stats.collect{it[1]})

    //
    // Create Gemini-compatible database files
    //

    if(params.gemini){

        CustomChannelOperators.joinOnKeys(
            annotation_output,
            VCF_EXTRACT_RELATE_SOMALIER.out.samples_tsv,
            ['id', 'family', 'family_count']
        )
        .dump(tag:'vcf2db_input', pretty:true)
        .set { vcf2db_input }

        VCF2DB(
            vcf2db_input
        )

        ch_versions = ch_versions.mix(VCF2DB.out.versions)

        VCF2DB.out.db
            .dump(tag:'vcf2db_output', pretty:true)

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

def get_family_id_from_ped(ped_file){

    // Check if there is a file
    if (ped_file.isEmpty()){
        return null
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

    // get family_id
    return (ped =~ /\n([^#]\w+)/)[0][1]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
