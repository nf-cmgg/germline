/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateAndConvertSamplesheet } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
ch_multiqc_logo     = params.multiqc_logo   ? file(params.multiqc_logo, checkIfExists: true)   : file("$projectDir/assets/CMGG_logo.png", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLE_PREPARATION            } from '../subworkflows/local/sample_preparation'
include { GERMLINE_VARIANT_CALLING      } from '../subworkflows/local/germline_variant_calling'
include { JOINT_GENOTYPING              } from '../subworkflows/local/joint_genotyping'
include { ANNOTATION                    } from '../subworkflows/local/annotation'
include { ADD_PED_HEADER                } from '../subworkflows/local/add_ped_header'
include { VCF_VALIDATE_SMALL_VARIANTS   } from '../subworkflows/local/vcf_validate_small_variants/main'

include { VCF_EXTRACT_RELATE_SOMALIER   } from '../subworkflows/nf-core/vcf_extract_relate_somalier/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX as FAIDX                                    } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_COMPOSESTRTABLEFILE as COMPOSESTRTABLEFILE           } from '../modules/nf-core/gatk4/composestrtablefile/main'
include { RTGTOOLS_FORMAT                                            } from '../modules/nf-core/rtgtools/format/main'
include { UNTAR                                                      } from '../modules/nf-core/untar/main'
include { TABIX_TABIX as TABIX_DBSNP                                 } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TRUTH                                 } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GVCF                                  } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_FINAL                                 } from '../modules/nf-core/tabix/tabix/main'
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

    //
    // Check input path parameters to see if they exist
    //

    def checkPathParamList = [
        params.fasta,
        params.fai,
        params.dict,
        params.strtablefile,
        params.dbsnp,
        params.dbsnp_tbi,
        params.somalier_sites,
        params.sdf,
        params.roi,
        params.vep_cache,
        params.vcfanno_config,
        params.vcfanno_lua,
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.spliceai_indel,
        params.spliceai_indel_tbi,
        params.spliceai_snv,
        params.spliceai_snv_tbi,
        params.mastermind,
        params.mastermind_tbi,
        params.maxentscan,
        params.eog,
        params.eog_tbi
    ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    //
    // Check the input samplesheet
    //

    if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { error('Input samplesheet not specified!') }


    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Importing and convert the input files passed through the parameters to channels
    //

    ch_fasta_ready        = Channel.fromPath(params.fasta).collect()
    ch_fai                = params.fai                 ? Channel.fromPath(params.fai).collect()                                        : null
    ch_dict               = params.dict                ? Channel.fromPath(params.dict).collect()                                       : null
    ch_strtablefile       = params.strtablefile        ? Channel.fromPath(params.strtablefile).collect()                               : null
    ch_sdf                = params.sdf                 ? Channel.fromPath(params.sdf).map {sdf -> [[id:'reference'], sdf]}.collect()   : null

    ch_default_roi        = params.roi                 ? Channel.fromPath(params.roi).collect()                : []

    ch_dbsnp_ready        = params.dbsnp               ? Channel.fromPath(params.dbsnp).collect()              : []
    ch_dbsnp_tbi          = params.dbsnp_tbi           ? Channel.fromPath(params.dbsnp_tbi).collect()          : []

    ch_somalier_sites     = params.somalier_sites      ? Channel.fromPath(params.somalier_sites).collect()     : []

    ch_vep_cache          = params.vep_cache           ? Channel.fromPath(params.vep_cache).collect()          : []

    ch_vcfanno_config     = params.vcfanno_config      ? Channel.fromPath(params.vcfanno_config).collect()     : []
    ch_vcfanno_lua        = params.vcfanno_lua         ? Channel.fromPath(params.vcfanno_lua).collect()        : []
    ch_vcfanno_resources  = params.vcfanno_resources   ? Channel.of(params.vcfanno_resources.split(",")).map({ file(it, checkIfExists:true) }).collect()   : []

    def val_callers       = params.callers.tokenize(",")
    def allowed_callers   = ["haplotypecaller", "deepvariant"]

    for (caller in val_callers) { if(!(caller in allowed_callers)) {error("The caller '${caller}' is not supported please specify a comma delimited list with on or more of the following callers: ${allowed_callers}".toString())}}

    //
    // Check for the presence of EnsemblVEP plugins that use extra files
    //

    ch_vep_extra_files = []

    if(params.annotate){
        // Check if all dbnsfp files are given
        if (params.dbnsfp && params.dbnsfp_tbi && params.vep_dbnsfp) {
            ch_vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
            ch_vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
        }
        else if (params.vep_dbnsfp) {
            exit 1, "Please specify '--vep_dbsnfp true', '--dbnsfp PATH/TO/DBNSFP/FILE' and '--dbnspf_tbi PATH/TO/DBNSFP/INDEX/FILE' to use the dbnsfp VEP plugin."
        }

        // Check if all spliceai files are given
        if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi && params.vep_spliceai) {
            ch_vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
            ch_vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
            ch_vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
            ch_vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        }
        else if (params.vep_spliceai) {
            exit 1, "Please specify '--vep_spliceai true', '--spliceai_snv PATH/TO/SPLICEAI/SNV/FILE', '--spliceai_snv_tbi PATH/TO/SPLICEAI/SNV/INDEX/FILE', '--spliceai_indel PATH/TO/SPLICEAI/INDEL/FILE' and '--spliceai_indel_tbi PATH/TO/SPLICEAI/INDEL/INDEX/FILE' to use the SpliceAI VEP plugin."
        }

        // Check if all mastermind files are given
        if (params.mastermind && params.mastermind_tbi && params.vep_mastermind) {
            ch_vep_extra_files.add(file(params.mastermind, checkIfExists: true))
            ch_vep_extra_files.add(file(params.mastermind_tbi, checkIfExists: true))
        }
        else if (params.vep_mastermind) {
            exit 1, "Please specify '--vep_mastermind true', '--mastermind PATH/TO/MASTERMIND/FILE' and '--mastermind_tbi PATH/TO/MASTERMIND/INDEX/FILE' to use the mastermind VEP plugin."
        }

        // Check if all maxentscan files are given
        if (params.maxentscan && params.vep_maxentscan) {
            ch_vep_extra_files.add(file(params.maxentscan, checkIfExists: true))
        }
        else if (params.vep_maxentscan) {
            exit 1, "Please specify '--vep_maxentscan true', '--maxentscan PATH/TO/MAXENTSCAN/' to use the MaxEntScan VEP plugin."
        }

        // Check if all EOG files are given
        if (params.eog && params.eog_tbi && params.vep_eog) {
            ch_vep_extra_files.add(file(params.eog, checkIfExists: true))
            ch_vep_extra_files.add(file(params.eog_tbi, checkIfExists: true))
        }
        else if (params.vep_eog) {
            exit 1, "Please specify '--vep_eog true', '--eog PATH/TO/EOG/FILE' and '--eog_tbi PATH/TO/EOG/INDEX/FILE' to use the EOG custom VEP plugin."
        }
    }

    //
    // Create the optional input files if they are not supplied
    //

    // DBSNP index
    if (ch_dbsnp_ready && !ch_dbsnp_tbi) {
        TABIX_DBSNP(
            ch_dbsnp_ready.map { [[id:'dbsnp'], it] }
        )
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)

        TABIX_DBSNP.out.tbi
            .map(
                { meta, tbi ->
                    [ tbi ]
                }
            )
            .collect()
            .set { ch_dbsnp_tbi_ready }
    } else if (ch_dbsnp_ready) {
        ch_dbsnp_tbi.set { ch_dbsnp_tbi_ready }
    } else {
        ch_dbsnp_tbi_ready = []
    }

    // Reference fasta index
    if (!ch_fai) {
        FAIDX(
            ch_fasta_ready.map({ fasta -> [ [id:"fasta_fai"], fasta ]})
        )
        ch_versions = ch_versions.mix(FAIDX.out.versions)

        FAIDX.out.fai
            .map({ meta, fai -> [ fai ]})
            .collect()
            .dump(tag:'fasta_fai', pretty:true)
            .set { ch_fai_ready }
    } else {
        ch_fai.set { ch_fai_ready }
    }

    // Reference sequence dictionary
    if (!ch_dict) {
        CREATESEQUENCEDICTIONARY(
            ch_fasta_ready
        )
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)

        CREATESEQUENCEDICTIONARY.out.dict
            .collect()
            .dump(tag:'dict', pretty:true)
            .set { ch_dict_ready }
    } else {
        ch_dict.set { ch_dict_ready }
    }

    // Reference STR table file
    if (params.dragstr && !ch_strtablefile) {
        COMPOSESTRTABLEFILE(
            ch_fasta_ready,
            ch_fai_ready,
            ch_dict_ready
        )
        ch_versions  = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions)

        COMPOSESTRTABLEFILE.out.str_table
            .collect()
            .dump(tag:'strtablefile', pretty:true)
            .set { ch_strtablefile_ready }
    } else if (params.dragstr) {
        ch_strtablefile.set { ch_strtablefile_ready }
    } else {
        ch_strtablefile_ready = []
    }

    // Reference validation SDF
    if (params.validate && !ch_sdf) {
        RTGTOOLS_FORMAT(
            ch_fasta_ready.map { fasta -> [[id:'reference'], fasta, [], []]}
        )
        ch_versions  = ch_versions.mix(RTGTOOLS_FORMAT.out.versions)

        RTGTOOLS_FORMAT.out.sdf
            .collect()
            .dump(tag:'sdf', pretty:true)
            .set { ch_sdf_ready }
    }
    else if (params.validate && params.sdf.endsWith(".tar.gz")) {
        UNTAR(
            ch_sdf
        )
        ch_versions = ch_versions.mix(UNTAR.out.versions)

        UNTAR.out.untar
            .dump(tag:'sdf', pretty:true)
            .set { ch_sdf_ready }
    } else if(params.validate) {
        ch_sdf.set { ch_sdf_ready }
    } else {
        ch_sdf_ready = [[],[]]
    }

    //
    // Read in samplesheet, validate and convert to a channel
    //

    Channel.validateAndConvertSamplesheet(
        file(params.input, checkIfExists:true),
        file("${projectDir}/assets/schema_input.json", checkIfExists:true)
    )
        .map { meta, cram, crai, gvcf, tbi, roi, ped, truth_vcf, truth_tbi, truth_bed ->
            // Infer the family ID from the PED file if no family ID was given.
            // If no PED is given, use the sample ID as family ID            
            new_meta = meta + [
                family: meta.family ?: ped ? get_family_id_from_ped(ped) : meta.sample, 
            ]
            [ new_meta, cram, crai, gvcf, tbi, roi, ped, truth_vcf, truth_tbi, truth_bed ]
        }
        .tap { ch_raw_inputs }
        .map { it[0] }
        .distinct() // Make sure the same sample isn't counted twice when given multiple times
        .map { meta ->
            meta.family
        }
        .reduce([:]) { counts, v ->
            // Count the unique samples in one family
            counts[v] = (counts[v] ?: 0) + 1
            counts
        }
        .combine(ch_raw_inputs)
        .multiMap { counts, meta, cram, crai, gvcf, tbi, roi, ped, truth_vcf, truth_tbi, truth_bed ->
            // Divide the input files into their corresponding channel
            new_meta_family = [
                id:             meta.family,
                family:         meta.family,
                family_count:   counts[meta.family] // Contains the amount of samples in the current family
            ]

            new_meta = meta + [
                family_count:   counts[meta.family], // Contains the amount of samples in the family from this sample
                type:           gvcf ? "gvcf" : "cram", // Whether a GVCF is given to this sample or not (aka skip variantcalling or not)
                analysis_type:  roi ? "WES" : "WGS" // Define if the data is WGS or WES
            ]

            truth_variants: [new_meta_family, truth_vcf, truth_tbi, truth_bed, meta.id] // Optional channel containing the truth VCF, its index and the optional BED file
            gvcf:           [new_meta, gvcf, tbi] // Optional channel containing the GVCFs and their optional indices
            cram:           [new_meta, cram, crai]  // Mandatory channel containing the CRAM files and their optional indices
            peds:           [new_meta_family, ped] // Optional channel containing the PED files 
            roi:            [new_meta, roi] // Optional channel containing the ROI BED files for WES samples
        }
        .set { ch_input }

    ch_input.gvcf.dump(tag:'input_gvcf', pretty:true)
    ch_input.roi.dump(tag:'input_roi', pretty:true)
    ch_input.truth_variants.dump(tag:'truth_variants', pretty:true)
    ch_input.cram.dump(tag:'input_crams', pretty:true)
    ch_input.peds.dump(tag:'input_peds', pretty:true)

    //
    // Create the GVCF index if it's missing
    //

    ch_input.gvcf
        .filter { it[0].type == "gvcf" } // Filter out samples which have no GVCF
        .branch { meta, gvcf, tbi ->
            no_tbi: !tbi
                return [ meta, gvcf ]
            tbi:    tbi
                return [ meta, gvcf, tbi ]
        }
        .set { ch_gvcf_branch }

    TABIX_GVCF(
        ch_gvcf_branch.no_tbi
    )
    ch_versions = ch_versions.mix(TABIX_GVCF.out.versions)

    // TODO find a clever solution for the possible double analysis of GVCFs (will only happen when two callers are used)
    ch_gvcf_branch.no_tbi
        .map { meta, vcf, tbi ->
            [ meta, vcf, tbi, val_callers ]
        }
        .transpose(by:3)
        .map { meta, vcf, tbi, caller ->
            new_meta = meta + [caller:caller]
            [ new_meta, vcf, tbi ]
        }
        .join(TABIX_GVCF.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_gvcf_branch.tbi)
        .set { ch_gvcfs_ready }

    //
    // Run sample preparation
    //

    SAMPLE_PREPARATION(
        ch_input.cram,
        ch_input.roi,
        ch_fasta_ready,
        ch_fai_ready,
        ch_default_roi
    )
    ch_versions = ch_versions.mix(SAMPLE_PREPARATION.out.versions)

    //
    // Take one PED file per family
    //

    ch_input.peds
        .groupTuple() // No size needed here because no process has been run with PED files before this
        .map { meta, peds ->
            // Find the first PED file and return that one for the family ([] if no PED is given for the family)
            [ meta, peds.find { it != [] } ?: [], val_callers ]
        }
        .transpose(by:2)
        .map { meta, ped, caller ->
            new_meta = meta + [caller:caller]
            [ new_meta, ped ]
        }
        .dump(tag:'peds', pretty:true)
        .set { ch_peds_ready }

    //
    // Perform the variant calling
    //

    GERMLINE_VARIANT_CALLING(
        SAMPLE_PREPARATION.out.ready_crams.filter { it[0].type == "cram" }, // Filter out files that already have a called GVCF
        SAMPLE_PREPARATION.out.ready_beds.filter { it[0].type == "cram" }, // Filter out files that already have a called GVCF
        ch_fasta_ready,
        ch_fai_ready,
        ch_dict_ready,
        ch_strtablefile_ready,
        ch_dbsnp_ready,
        ch_dbsnp_tbi_ready
    )
    ch_versions = ch_versions.mix(GERMLINE_VARIANT_CALLING.out.versions)
    ch_reports  = ch_reports.mix(GERMLINE_VARIANT_CALLING.out.reports)

    GERMLINE_VARIANT_CALLING.out.gvcfs
        .mix(ch_gvcfs_ready)
        .dump(tag:'variantcalling_output', pretty:true)
        .set { ch_variantcalling_output }

    //
    // Joint-genotyping of the families
    //

    JOINT_GENOTYPING(
        ch_variantcalling_output,
        SAMPLE_PREPARATION.out.ready_beds,
        ch_fasta_ready,
        ch_fai_ready,
        ch_dict_ready,
        ch_dbsnp_ready,
        ch_dbsnp_tbi_ready
    )
    ch_versions = ch_versions.mix(JOINT_GENOTYPING.out.versions)

    JOINT_GENOTYPING.out.genotyped_vcfs
        .dump(tag:'joint_genotyping_output', pretty:true)
        .set { ch_joint_genotyping_output }

    //
    // Filter the variants
    //

    if (params.filter) {
        FILTER_SNPS(
            ch_joint_genotyping_output
        )
        ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)

        FILTER_INDELS(
            FILTER_SNPS.out.vcf
        )
        ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

        FILTER_INDELS.out.vcf.set { ch_filter_output }
    } else {
        ch_joint_genotyping_output.set { ch_filter_output }
    }

    ch_filter_output.dump(tag:'filter_output', pretty: true)

    //
    // Run relation tests with somalier
    //

    VCF_EXTRACT_RELATE_SOMALIER(
        ch_filter_output
            .filter { meta, vcf ->
                // Filter out the families that only have one individual
                meta.family_count > 1
            }
            .map { it + [[], 1] },
        ch_fasta_ready,
        ch_fai_ready,
        ch_somalier_sites,
        ch_peds_ready
            .filter { meta, ped ->
                // Filter out the families that only have one individual
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
        ch_filter_output
            .branch { meta, vcf ->
                // Only add ped headers to VCFs with more than one individual
                single: meta.family_count == 1
                multiple: meta.family_count > 1
            }
            .set { ch_ped_header_branch }

        ADD_PED_HEADER(
            ch_ped_header_branch.multiple,
            VCF_EXTRACT_RELATE_SOMALIER.out.samples_tsv
        )
        ch_versions = ch_versions.mix(ADD_PED_HEADER.out.versions)

        ADD_PED_HEADER.out.ped_vcfs
            .mix(ch_ped_header_branch.single)
            .dump(tag:'ped_vcfs', pretty:true)
            .set { ch_ped_vcfs }
    } else {
        ch_filter_output.set { ch_ped_vcfs }
    }

    //
    // Annotation of the variants and creation of Gemini-compatible database files
    //

    if (params.annotate) {
        ANNOTATION(
            ch_ped_vcfs,
            ch_fasta_ready,
            ch_fai_ready,
            ch_vep_cache,
            ch_vep_extra_files,
            ch_vcfanno_config,
            ch_vcfanno_lua,
            ch_vcfanno_resources
        )
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)
        ch_reports  = ch_reports.mix(ANNOTATION.out.reports)

        ANNOTATION.out.annotated_vcfs.set { ch_annotation_output }
    } else {
        ch_ped_vcfs.set { ch_annotation_output }
    }

    ch_annotation_output.dump(tag:'annotation_output', pretty:true)

    //
    // Tabix the resulting VCF
    //

    TABIX_FINAL(
        ch_annotation_output
    )
    ch_versions = ch_versions.mix(TABIX_FINAL.out.versions.first())

    ch_annotation_output
        .join(TABIX_FINAL.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_final_vcfs }

    //
    // Validate the found variants
    //

    if (params.validate){

        ch_input.truth_variants
            .map { meta, vcf, tbi, bed, sample ->
                [ meta, vcf, tbi, bed, sample, val_callers ]
            }
            .transpose(by:5)
            .map { meta, vcf, tbi, bed, sample, caller ->
                new_meta = meta + [caller:caller]
                [ new_meta, vcf, tbi, bed, sample ]
            }
            .set { ch_truth_variants_ready }

        ch_final_vcfs
            .combine(ch_truth_variants_ready, by: 0)
            .map { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed, sample ->
                new_meta = meta + [sample:sample]
                [ new_meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ]
            }
            .filter { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ->
                // Filter out all samples that have no truth VCF
                truth_vcf != []
            }
            .branch { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ->
                tbi: truth_tbi != []
                no_tbi: truth_tbi == []
            }
            .set { ch_validation_branch }

        // Create truth VCF indices if none were given
        TABIX_TRUTH(
            ch_validation_branch.no_tbi.map { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed -> 
                [ meta, truth_vcf ]
            }
        )
        ch_versions = ch_versions.mix(TABIX_TRUTH.out.versions)

        ch_validation_branch.no_tbi
            .join(TABIX_TRUTH.out.tbi, failOnDuplicate: true, failOnMismatch: true) 
            .map { meta, vcf, tbi, truth_vcf, empty_tbi, truth_bed, truth_tbi ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ]
            }
            .mix(ch_validation_branch.tbi)
            .multiMap { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ->
                vcfs: [meta, vcf, tbi, truth_vcf, truth_tbi]
                bed:  [meta, truth_bed, []]
            }
            .set { ch_validation_input }

        VCF_VALIDATE_SMALL_VARIANTS(
            ch_validation_input.vcfs,
            ch_validation_input.bed,
            ch_fasta_ready.map { [[], it] },
            ch_fai_ready.map { [[], it] },
            ch_sdf_ready.collect(),
            [[],[]],
            [[],[]],
            [[],[]],
            "vcfeval" //Only VCFeval for now, awaiting the conda fix for happy (https://github.com/bioconda/bioconda-recipes/pull/39267)
        )
        ch_versions = ch_versions.mix(VCF_VALIDATE_SMALL_VARIANTS.out.versions)
    }

    //
    // Perform QC on the final VCFs
    //

    BCFTOOLS_STATS_FAMILY(
        ch_final_vcfs,
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
            ch_final_vcfs.map { meta, vcf, tbi -> [ meta, vcf ]},
            VCF_EXTRACT_RELATE_SOMALIER.out.samples_tsv,
            ['id', 'family', 'family_count']
        )
        .dump(tag:'vcf2db_input', pretty:true)
        .set { ch_vcf2db_input }

        VCF2DB(
            ch_vcf2db_input
        )
        ch_versions = ch_versions.mix(VCF2DB.out.versions.first())

        VCF2DB.out.db.dump(tag:'vcf2db_output', pretty:true)
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
        ch_multiqc_logo
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
