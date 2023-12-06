/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet } from 'plugin/nf-validation'

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
    if (!params.genome) { exit 1, "A genome should be supplied for annotation (use --genome)"}

    // Check if the VEP versions were given
    if (!params.vep_version) { exit 1, "A VEP version should be supplied for annotation (use --vep_version)"}
    if (!params.vep_cache_version) { exit 1, "A VEP cache version should be supplied for annotation (use --vep_cache_version)"}

    // Check if a species is entered
    if (!params.species) { exit 1, "A species should be supplied for annotation (use --species)"}

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

include { CRAM_PREPARE_SAMTOOLS_BEDTOOLS    } from '../subworkflows/local/cram_prepare_samtools_bedtools/main'
include { INPUT_SPLIT_BEDTOOLS              } from '../subworkflows/local/input_split_bedtools/main'
include { CRAM_CALL_GENOTYPE_GATK4          } from '../subworkflows/local/cram_call_genotype_gatk4/main'
include { CRAM_CALL_VARDICTJAVA             } from '../subworkflows/local/cram_call_vardictjava/main'
include { VCF_EXTRACT_RELATE_SOMALIER       } from '../subworkflows/local/vcf_extract_relate_somalier/main'
include { VCF_PED_RTGTOOLS                  } from '../subworkflows/local/vcf_ped_rtgtools/main'
include { VCF_ANNOTATION                    } from '../subworkflows/local/vcf_annotation/main'
include { VCF_VALIDATE_SMALL_VARIANTS       } from '../subworkflows/local/vcf_validate_small_variants/main'

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
include { ENSEMBLVEP_DOWNLOAD                                        } from '../modules/nf-core/ensemblvep/download/main'
include { BCFTOOLS_STATS                                             } from '../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_NORM                                              } from '../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_DECOMPOSE                             } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_NORMALIZE                             } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DBSNP                                 } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GVCF                                  } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TRUTH                                 } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_FINAL                                 } from '../modules/nf-core/tabix/tabix/main'
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
        params.eog,
        params.eog_tbi
    ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    //
    // Check the input samplesheet
    //

    if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { error('Input samplesheet not specified!') }

    callers = params.callers.tokenize(",")
    for(caller in callers) {
        if(!(caller in GlobalVariables.availableCallers)) { error("\"${caller}\" is not a supported callers please use one or more of these instead: ${GlobalVariables.availableCallers}")}
    }

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Importing and convert the input files passed through the parameters to channels
    //

    ch_fasta_ready        = Channel.fromPath(params.fasta).map{ [[id:"reference"], it]}.collect()
    ch_fai                = params.fai                 ? Channel.fromPath(params.fai).map{ [[id:"reference"], it]}.collect()                                        : null
    ch_dict               = params.dict                ? Channel.fromPath(params.dict).map{ [[id:"reference"], it]}.collect()                                       : null
    ch_strtablefile       = params.strtablefile        ? Channel.fromPath(params.strtablefile).map{ [[id:"reference"], it]}.collect()                               : null
    ch_sdf                = params.sdf                 ? Channel.fromPath(params.sdf).map {sdf -> [[id:'reference'], sdf]}.collect()   : null

    ch_default_roi        = params.roi                 ? Channel.fromPath(params.roi).collect()                : []

    ch_dbsnp_ready        = params.dbsnp               ? Channel.fromPath(params.dbsnp).collect()              : Channel.value([])
    ch_dbsnp_tbi          = params.dbsnp_tbi           ? Channel.fromPath(params.dbsnp_tbi).collect()          : Channel.value([])

    ch_somalier_sites     = params.somalier_sites      ? Channel.fromPath(params.somalier_sites).collect()     : []

    ch_vep_cache          = params.vep_cache           ? Channel.fromPath(params.vep_cache).collect()          : []

    ch_vcfanno_config     = params.vcfanno_config      ? Channel.fromPath(params.vcfanno_config).collect()     : []
    ch_vcfanno_lua        = params.vcfanno_lua         ? Channel.fromPath(params.vcfanno_lua).collect()        : []
    ch_vcfanno_resources  = params.vcfanno_resources   ? Channel.of(params.vcfanno_resources.split(",")).map({ file(it, checkIfExists:true) }).collect()   : []

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
            .map{ meta, tbi ->
                [ tbi ]
            }
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
            ch_fasta_ready
        )
        ch_versions = ch_versions.mix(FAIDX.out.versions)

        FAIDX.out.fai
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
            ch_fasta_ready.map { meta, fasta -> [meta, fasta, [], []]}
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

    if (!ch_vep_cache && params.annotate) {
        ENSEMBLVEP_DOWNLOAD(
            Channel.of([[id:"vep_cache"], params.genome == "hg38" ? "GRCh38" : params.genome, params.species, params.vep_cache_version]).collect()
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

        ch_vep_cache_ready = ENSEMBLVEP_DOWNLOAD.out.cache.map{it[1]}.collect()
    } else {
        ch_vep_cache_ready = ch_vep_cache
    }

    //
    // Read in samplesheet, validate and convert to a channel
    //

    // Output the samplesheet
    file(params.input).copyTo("${params.outdir}/samplesheet.csv")

    Channel.fromSamplesheet("input", immutable_meta: false)
        .map { meta, cram, crai, gvcf, tbi, roi, ped, truth_vcf, truth_tbi, truth_bed ->
            // Infer the family ID from the PED file if no family ID was given.
            // If no PED is given, use the sample ID as family ID            
            def new_meta = meta + [
                family: meta.family ?: ped ? get_family_id_from_ped(ped) : meta.sample, 
            ]
            [ new_meta, cram, crai, gvcf, tbi, roi, ped, truth_vcf, truth_tbi, truth_bed ]
        }
        .tap { ch_raw_inputs }
        .map { [ "id":it[0].id, "family":it[0].family ] }
        .reduce([:]) { families, v ->
            // Count the unique samples in one family
            families[v.family] = families[v.family] ? families[v.family] + [v.id] : [v.id]
            families[v.family] = families[v.family].unique()
            families
        }
        .combine(ch_raw_inputs)
        .multiMap { families, meta, cram, crai, gvcf, tbi, roi, ped, truth_vcf, truth_tbi, truth_bed ->
            // Divide the input files into their corresponding channel
            def new_meta = meta + [
                family_count:   families[meta.family].size(), // Contains the amount of samples in the family from this sample
                type: gvcf && cram ? "gvcf_cram" : gvcf ? "gvcf" : "cram" // Define the type of input data
            ]

            def new_meta_ped = meta - meta.subMap(["type", "family_count", "vardict_min_af"])

            def new_meta_validation = [
                id: meta.id,
                sample: meta.sample,
                family: meta.family
            ]

            truth_variants: [new_meta_validation, truth_vcf, truth_tbi, truth_bed] // Optional channel containing the truth VCF, its index and the optional BED file
            gvcf:           [new_meta, gvcf, tbi] // Optional channel containing the GVCFs and their optional indices
            cram:           [new_meta, cram, crai]  // Mandatory channel containing the CRAM files and their optional indices
            peds:           [new_meta_ped, ped] // Optional channel containing the PED files 
            roi:            [new_meta, roi] // Optional channel containing the ROI BED files for WES samples
            family_samples: [meta.family, families[meta.family]] // A channel containing the samples per family
        }
        .set { ch_input }

    ch_family_samples = ch_input.family_samples.distinct()

    //
    // Create the GVCF index if it's missing
    //

    ch_input.gvcf
        .filter { it[0].type == "gvcf" || it[0].type == "gvcf_cram" } // Filter out samples that have no GVCF
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

    ch_gvcf_branch.no_tbi
        .join(TABIX_GVCF.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_gvcf_branch.tbi)
        .set { ch_gvcfs_ready }

    //
    // Run sample preparation
    //

    CRAM_PREPARE_SAMTOOLS_BEDTOOLS(
        ch_input.cram.filter { it[0].type == "cram" || (it[0].type == "gvcf_cram" && callers - GlobalVariables.gvcfCallers) }, // Filter out files that already have a called GVCF when only GVCF callers are used
        ch_input.roi.filter { it[0].type == "cram" || (it[0].type == "gvcf_cram" && callers - GlobalVariables.gvcfCallers) }, // Filter out files that already have a called GVCF when only GVCF callers are used
        ch_fasta_ready,
        ch_fai_ready,
        ch_default_roi
    )
    ch_versions = ch_versions.mix(CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.versions)

    //
    // Split the BED files
    //

    INPUT_SPLIT_BEDTOOLS(
        CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_beds.map { it + [params.scatter_count] },
        CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams
    )
    ch_versions = ch_versions.mix(INPUT_SPLIT_BEDTOOLS.out.versions)

    ch_calls = Channel.empty()

    if("haplotypecaller" in callers) {
            
        //
        // Call variants with GATK4 HaplotypeCaller
        //

        CRAM_CALL_GENOTYPE_GATK4(
            INPUT_SPLIT_BEDTOOLS.out.split.filter { it[0].type == "cram" }, // Filter out the entries that already have a GVCF
            ch_gvcfs_ready,
            ch_fasta_ready,
            ch_fai_ready,
            ch_dict_ready,
            ch_strtablefile_ready,
            ch_dbsnp_ready,
            ch_dbsnp_tbi_ready
        )
        ch_versions = ch_versions.mix(CRAM_CALL_GENOTYPE_GATK4.out.versions)
        ch_reports  = ch_reports.mix(CRAM_CALL_GENOTYPE_GATK4.out.reports)

        ch_calls = ch_calls.mix(CRAM_CALL_GENOTYPE_GATK4.out.vcfs)

    }

    if("vardict" in callers) {
            
        //
        // Call variants with VarDict
        //

        CRAM_CALL_VARDICTJAVA(
            CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams,
            INPUT_SPLIT_BEDTOOLS.out.split,
            ch_fasta_ready,
            ch_fai_ready
        )
        ch_versions = ch_versions.mix(CRAM_CALL_VARDICTJAVA.out.versions)

        ch_calls = ch_calls.mix(CRAM_CALL_VARDICTJAVA.out.vcfs)
    
    }

    ch_calls
        .map { meta, vcf, tbi ->
            def new_meta = meta - meta.subMap(["type", "vardict_min_af"])
            [ new_meta, vcf, tbi ]
        }
        .set { ch_called_variants }

    BCFTOOLS_STATS(
        ch_called_variants,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
    ch_reports = ch_reports.mix(BCFTOOLS_STATS.out.stats.collect { it[1] })

    if(params.normalize) {
        BCFTOOLS_NORM(
            ch_called_variants,
            ch_fasta_ready,
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

        TABIX_NORMALIZE(
            BCFTOOLS_NORM.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_NORMALIZE.out.versions.first())

        BCFTOOLS_NORM.out.vcf
            .join(TABIX_NORMALIZE.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { ch_normalized_variants }
    } else {
        ch_called_variants.set { ch_normalized_variants }
    }

    if(!params.only_merge && !params.only_call) {

        //
        // Preprocess the PED channel
        //

        ch_input.peds
            .map { meta, ped ->
                [ meta.family, ped ]
            }
            .groupTuple() // No size needed here because no process has been run with PED files before this
            .map { meta, peds ->
                // Find the first PED file and return that one for the family ([] if no PED is given for the family)
                [ meta, peds.find { it != [] } ?: [] ]
            }
            .combine(ch_normalized_variants.map { meta, vcf, tbi -> [ meta.family, meta, vcf, tbi ]}, by:0)
            .map { family, ped, meta, vcf, tbi ->
                [ meta, ped ]
            }
            .set { ch_somalier_input }

        //
        // Run relation tests with somalier
        //

        VCF_EXTRACT_RELATE_SOMALIER(
            ch_normalized_variants,
            ch_fasta_ready.map { it[1] },
            ch_fai_ready.map { it[1] },
            ch_somalier_sites,
            ch_somalier_input
        )
        ch_versions = ch_versions.mix(VCF_EXTRACT_RELATE_SOMALIER.out.versions)

        //
        // Add PED headers to the VCFs
        //

        if(params.add_ped){

            VCF_PED_RTGTOOLS(
                ch_normalized_variants,
                VCF_EXTRACT_RELATE_SOMALIER.out.peds
            )
            ch_versions = ch_versions.mix(VCF_PED_RTGTOOLS.out.versions)

            VCF_PED_RTGTOOLS.out.ped_vcfs
                .set { ch_ped_vcfs }
        } else {
            ch_normalized_variants
                .map { meta, vcf, tbi=[] ->
                    [ meta, vcf ]
                }
                .set { ch_ped_vcfs }
        }

        //
        // Annotation of the variants and creation of Gemini-compatible database files
        //

        if (params.annotate) {
            VCF_ANNOTATION(
                ch_ped_vcfs,
                ch_fasta_ready,
                ch_fai_ready,
                ch_vep_cache_ready,
                ch_vep_extra_files,
                ch_vcfanno_config,
                ch_vcfanno_lua,
                ch_vcfanno_resources
            )
            ch_versions = ch_versions.mix(VCF_ANNOTATION.out.versions)
            ch_reports  = ch_reports.mix(VCF_ANNOTATION.out.reports)

            VCF_ANNOTATION.out.annotated_vcfs.set { ch_annotation_output }
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
                .groupTuple() // No size needed here since it's being run before any process
                .map { meta, vcf, tbi, bed ->
                    // Get only one VCF for samples that were given multiple times
                    one_vcf = vcf.find { it != [] } ?: []
                    one_tbi = tbi.find { it != [] } ?: []
                    one_bed = bed.find { it != [] } ?: []
                    [ meta, one_vcf, one_tbi, one_bed ]
                }
                .branch { meta, vcf, tbi, bed ->
                    no_vcf: !vcf
                    tbi: tbi
                    no_tbi: !tbi
                }
                .set { ch_truths_input }

            // Create truth VCF indices if none were given
            TABIX_TRUTH(
                ch_truths_input.no_tbi.map { meta, vcf, tbi, bed -> 
                    [ meta, vcf ]
                }
            )
            ch_versions = ch_versions.mix(TABIX_TRUTH.out.versions.first())         

            ch_truths_input.no_tbi
                .join(TABIX_TRUTH.out.tbi, failOnDuplicate:true, failOnMismatch:true)
                .map { meta, vcf, empty, bed, tbi ->
                    [ meta, vcf, tbi, bed ]
                }
                .mix(ch_truths_input.tbi)
                .mix(ch_truths_input.no_vcf)
                .combine(callers)
                .map { meta, vcf, tbi, bed, caller ->
                    def new_meta = meta + [caller: caller]
                    [ new_meta, vcf, tbi, bed ]
                }
                .set { ch_truths }

            ch_final_vcfs
                .map { meta, vcf, tbi ->
                    def new_meta = meta - meta.subMap("family_count")
                    [ meta.family, new_meta, vcf, tbi ]
                }
                .combine(ch_family_samples, by:0)
                .map { family, meta, vcf, tbi, samples ->
                    def sample = meta.sample ? [meta.sample] : samples
                    [ meta, vcf, tbi, sample ]
                }
                .transpose(by: 3)
                .map { meta, vcf, tbi, sample ->
                    def new_meta = [
                        id: sample,
                        sample: sample,
                        family: meta.family,
                        caller: meta.caller
                    ]
                    [ new_meta, vcf, tbi ]
                }
                .combine(ch_truths, by:0)
                .filter { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ->
                    // Filter out all samples that have no truth VCF
                    truth_vcf != []
                }
                .multiMap { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ->
                    vcfs: [meta, vcf, tbi, truth_vcf, truth_tbi]
                    bed:  [meta, truth_bed, []]
                }
                .set { ch_validation_input }

            VCF_VALIDATE_SMALL_VARIANTS(
                ch_validation_input.vcfs,
                ch_validation_input.bed,
                ch_fasta_ready,
                ch_fai_ready,
                ch_sdf_ready.collect(),
                [[],[]],
                [[],[]],
                [[],[]],
                "vcfeval" //Only VCFeval for now, awaiting the conda fix for happy (https://github.com/bioconda/bioconda-recipes/pull/39267)
            )
            ch_versions = ch_versions.mix(VCF_VALIDATE_SMALL_VARIANTS.out.versions)
        }

        //
        // Create Gemini-compatible database files
        //

        if(params.gemini){
            CustomChannelOperators.joinOnKeys(
                ch_final_vcfs.map { meta, vcf, tbi -> [ meta, vcf ]},
                VCF_EXTRACT_RELATE_SOMALIER.out.peds,
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
    NfcoreTemplate.dump_parameters(workflow, params)
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
