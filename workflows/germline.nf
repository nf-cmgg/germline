/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_cmgg_germline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CRAM_PREPARE_SAMTOOLS_BEDTOOLS    } from '../subworkflows/local/cram_prepare_samtools_bedtools/main'
include { INPUT_SPLIT_BEDTOOLS              } from '../subworkflows/local/input_split_bedtools/main'
include { CRAM_CALL_GATK4                   } from '../subworkflows/local/cram_call_gatk4/main'
include { GVCF_JOINT_GENOTYPE_GATK4         } from '../subworkflows/local/gvcf_joint_genotype_gatk4/main'
include { BAM_CALL_ELPREP                   } from '../subworkflows/local/bam_call_elprep/main'
include { BAM_CALL_VARDICTJAVA              } from '../subworkflows/local/bam_call_vardictjava/main'
include { VCF_EXTRACT_RELATE_SOMALIER       } from '../subworkflows/local/vcf_extract_relate_somalier/main'
include { VCF_PED_RTGTOOLS                  } from '../subworkflows/local/vcf_ped_rtgtools/main'
include { VCF_ANNOTATION                    } from '../subworkflows/local/vcf_annotation/main'
include { VCF_VALIDATE_SMALL_VARIANTS       } from '../subworkflows/local/vcf_validate_small_variants/main'
include { VCF_UPD_UPDIO                     } from '../subworkflows/local/vcf_upd_updio/main'
include { VCF_ROH_AUTOMAP                   } from '../subworkflows/local/vcf_roh_automap/main'
include { VCF_FILTER_BCFTOOLS               } from '../subworkflows/local/vcf_filter_bcftools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX as FAIDX                                    } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { ELPREP_FASTATOELFASTA                                      } from '../modules/nf-core/elprep/fastatoelfasta/main'
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
include { MULTIQC                                                    } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The main workflow
workflow GERMLINE {

    take:
    // Input channels
    ch_samplesheet              // queue channel: The input channel

    // File inputs
    fasta                       // string: path to the reference fasta
    fai                         // string: path to the index of the reference fasta
    dict                        // string: path to the sequence dictionary file
    elfasta                     // string: path to the elfasta reference file
    strtablefile                // string: path to the strtable file
    sdf                         // string: path to the SDF directory
    dbsnp                       // string: path to the DBSNP VCF file
    dbsnp_tbi                   // string: path to the index of the DBSNP VCF file
    vep_cache                   // string: path to the VEP cache
    dbnsfp                      // string: path to the DBNSFP file
    dbnsfp_tbi                  // string: path to the index of the DBNSFP file
    spliceai_indel              // string: path to the SpliceAI indels file
    spliceai_indel_tbi          // string: path to the index of the SpliceAI indels file
    spliceai_snv                // string: path to the SpliceAI SNV file
    spliceai_snv_tbi            // string: path to the index of the SpliceAI SNV file
    mastermind                  // string: path to the Mastermind file
    mastermind_tbi              // string: path to the index of the Mastermind file
    eog                         // string: path to the EOG file
    eog_tbi                     // string: path to the index of the EOG file
    alphamissense               // string: path to the Alphamissense file
    alphamissense_tbi           // string: path to the index of the Alphamissense file
    vcfanno_resources           // string: semicolon-separated paths/globs to the VCFanno resources
    vcfanno_config              // string: path to VCFanno config TOML
    multiqc_config              // string: path to the multiqc config file
    multiqc_logo                // string: path to the multiqc logo
    multiqc_methods_description // string: path to the multiqc methods description
    roi                         // string: path to the default ROI BED file
    somalier_sites              // string: path to the Somalier sites file
    vcfanno_lua                 // string: path to the VCFanno Lua script
    updio_common_cnvs           // string: path to the file containing common UPDio CNVs
    automap_repeats             // string: path to the Automap repeats file
    automap_panel               // string: path to the Automap panel file
    outdir                      // string: path to the output directory
    pedFiles                    // map:    a map that has the family ID as key and a PED file as value
    elsites                     // string: path to the elsites file for elprep

    // Boolean inputs
    dragstr                     // boolean: create a dragstr model and use it for haplotypecaller
    annotate                    // boolean: perform annotation
    vcfanno                     // boolean: use vcfanno annotations
    only_call                   // boolean: only perform variant calling
    only_merge                  // boolean: run the pipeline until after the family merge
    filter                      // boolean: filter the VCFs
    normalize                   // boolean: perform normalization
    add_ped                     // boolean: add ped headers to each VCF
    gemini                      // boolean: convert the VCF to a gemini database
    validate                    // boolean: validate the pipeline output when a truth set has been given
    updio                       // boolean: run UPDio on the final VCFs
    automap                     // boolean: run Automap on the final VCFs
    vep_dbnsfp                  // boolean: use the DBNSFP VEP plugin
    vep_spliceai                // boolean: use the SpliceAI VEP plugin
    vep_mastermind              // boolean: use the Mastermind VEP plugin
    vep_eog                     // boolean: use the EOG VEP plugin
    vep_alphamissense           // boolean: use the AlphaMissense VEP plugin

    // Value inputs
    genome                      // string:  the genome used by the pipeline run
    species                     // string:  the species used by the pipeline run
    vep_cache_version           // integer: the vep cache version to be used
    vep_chunk_size              // integer: the chunk size to split each VCF file into for VEP
    scatter_count               // integer: the amount of scattering performed on each file
    callers                     // list:    the callers to use


    main:
    def ch_versions      = Channel.empty()
    def ch_reports       = Channel.empty()
    def ch_multiqc_files = Channel.empty()

    //
    // Importing and convert the input files passed through the parameters to channels
    //

    def ch_fasta_ready        = Channel.fromPath(fasta).map{ fasta_file -> [[id:"reference"], fasta_file] }.collect()
    def ch_fai                = fai                 ? Channel.fromPath(fai).map{ fai_file -> [[id:"reference"], fai_file] }.collect() : null
    def ch_dict               = dict                ? Channel.fromPath(dict).map{ dict_file -> [[id:"reference"], dict_file] }.collect() : null
    def ch_elfasta            = elfasta             ? Channel.fromPath(elfasta).map { elfasta_file -> [[id:"reference"], elfasta_file]}.collect() : null
    def ch_strtablefile       = strtablefile        ? Channel.fromPath(strtablefile).map{ str_file -> [[id:"reference"], str_file] }.collect() : null
    def ch_sdf                = sdf                 ? Channel.fromPath(sdf).map { sdf_file -> [[id:'reference'], sdf_file] }.collect() : null

    def ch_default_roi        = roi                 ? Channel.fromPath(roi).collect()                : []

    def ch_dbsnp_ready        = dbsnp               ? Channel.fromPath(dbsnp).collect { dbsnp_file -> [[id:"dbsnp"], dbsnp_file] } : [[],[]]
    def ch_dbsnp_tbi          = dbsnp_tbi           ? Channel.fromPath(dbsnp_tbi).collect { dbsnp_file -> [[id:"dbsnp"], dbsnp_file] } : [[],[]]

    def ch_somalier_sites     = somalier_sites      ? Channel.fromPath(somalier_sites).collect { sites_file -> [[id:"somalier_sites"], sites_file] } : [[],[]]

    def ch_vep_cache          = vep_cache           ? Channel.fromPath(vep_cache).collect()          : []

    def ch_vcfanno_config     = vcfanno_config      ? Channel.fromPath(vcfanno_config).collect()     : []
    def ch_vcfanno_lua        = vcfanno_lua         ? Channel.fromPath(vcfanno_lua).collect()        : []
    def ch_vcfanno_resources  = vcfanno_resources   ? Channel.of(vcfanno_resources.split(";")).collect{ res -> file(res, checkIfExists:true) } : []

    def ch_updio_common_cnvs  = updio_common_cnvs   ? Channel.fromPath(updio_common_cnvs).map{ common_cnvs -> [[id:'updio_cnv'], common_cnvs] } : [[],[]]

    def ch_automap_repeats    = automap_repeats     ? Channel.fromPath(automap_repeats).map{ repeats ->  [[id:"repeats"], repeats] }.collect() : []
    def ch_automap_panel      = automap_panel       ? Channel.fromPath(automap_panel).map{ panel -> [[id:"automap_panel"], panel] }.collect() : [[],[]]

    def ch_elsites            = elsites             ? Channel.fromPath(elsites).map{ elsites_file -> [[id:'elsites'], elsites_file] }.collect() : [[],[]]

    //
    // Check for the presence of EnsemblVEP plugins that use extra files
    //

    def ch_vep_extra_files = []

    if(annotate){
        // Check if all dbnsfp files are given
        if (dbnsfp && dbnsfp_tbi && vep_dbnsfp) {
            ch_vep_extra_files.add(file(dbnsfp, checkIfExists: true))
            ch_vep_extra_files.add(file(dbnsfp_tbi, checkIfExists: true))
        }
        else if (vep_dbnsfp) {
            error("Please specify '--vep_dbsnfp true', '--dbnsfp PATH/TO/DBNSFP/FILE' and '--dbnspf_tbi PATH/TO/DBNSFP/INDEX/FILE' to use the dbnsfp VEP plugin.")
        }

        // Check if all spliceai files are given
        if (spliceai_snv && spliceai_snv_tbi && spliceai_indel && spliceai_indel_tbi && vep_spliceai) {
            ch_vep_extra_files.add(file(spliceai_snv, checkIfExists: true))
            ch_vep_extra_files.add(file(spliceai_snv_tbi, checkIfExists: true))
            ch_vep_extra_files.add(file(spliceai_indel, checkIfExists: true))
            ch_vep_extra_files.add(file(spliceai_indel_tbi, checkIfExists: true))
        }
        else if (vep_spliceai) {
            error("Please specify '--vep_spliceai true', '--spliceai_snv PATH/TO/SPLICEAI/SNV/FILE', '--spliceai_snv_tbi PATH/TO/SPLICEAI/SNV/INDEX/FILE', '--spliceai_indel PATH/TO/SPLICEAI/INDEL/FILE' and '--spliceai_indel_tbi PATH/TO/SPLICEAI/INDEL/INDEX/FILE' to use the SpliceAI VEP plugin.")
        }

        // Check if all mastermind files are given
        if (mastermind && mastermind_tbi && vep_mastermind) {
            ch_vep_extra_files.add(file(mastermind, checkIfExists: true))
            ch_vep_extra_files.add(file(mastermind_tbi, checkIfExists: true))
        }
        else if (vep_mastermind) {
            error("Please specify '--vep_mastermind true', '--mastermind PATH/TO/MASTERMIND/FILE' and '--mastermind_tbi PATH/TO/MASTERMIND/INDEX/FILE' to use the mastermind VEP plugin.")
        }

        // Check if all EOG files are given
        if (eog && eog_tbi && vep_eog) {
            ch_vep_extra_files.add(file(eog, checkIfExists: true))
            ch_vep_extra_files.add(file(eog_tbi, checkIfExists: true))
        }
        else if (vep_eog) {
            error("Please specify '--vep_eog true', '--eog PATH/TO/EOG/FILE' and '--eog_tbi PATH/TO/EOG/INDEX/FILE' to use the EOG custom VEP plugin.")
        }

        // Check if all AlphaMissense files are given
        if (alphamissense && alphamissense_tbi && vep_alphamissense) {
            ch_vep_extra_files.add(file(alphamissense, checkIfExists: true))
            ch_vep_extra_files.add(file(alphamissense_tbi, checkIfExists: true))
        }
        else if (vep_alphamissense) {
            error("Please specify '--vep_alphamissense true', '--alphamissense PATH/TO/ALPHAMISSENSE/FILE' and '--alphamissense_tbi PATH/TO/ALPHAMISSENSE/INDEX/FILE' to use the AlphaMissense VEP plugin.")
        }
    }

    //
    // Create the optional input files if they are not supplied
    //

    // DBSNP index
    def ch_dbsnp_tbi_ready = Channel.empty()
    if (ch_dbsnp_ready != [[],[]] && ch_dbsnp_tbi == [[],[]]) {
        TABIX_DBSNP(
            ch_dbsnp_ready
        )
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)

        ch_dbsnp_tbi_ready = TABIX_DBSNP.out.tbi.collect()
    } else {
        ch_dbsnp_tbi_ready = ch_dbsnp_tbi
    }

    // Reference fasta index
    def ch_fai_ready = Channel.empty()
    if (!ch_fai) {
        FAIDX(
            ch_fasta_ready,
            [[],[]]
        )
        ch_versions = ch_versions.mix(FAIDX.out.versions)

        ch_fai_ready = FAIDX.out.fai
            .collect()
    } else {
        ch_fai_ready = ch_fai
    }

    // Reference sequence dictionary
    def ch_dict_ready = Channel.empty()
    if (!ch_dict) {
        CREATESEQUENCEDICTIONARY(
            ch_fasta_ready
        )
        ch_versions = ch_versions.mix(CREATESEQUENCEDICTIONARY.out.versions)

        ch_dict_ready = CREATESEQUENCEDICTIONARY.out.dict
            .collect()
    } else {
        ch_dict_ready = ch_dict
    }

    def ch_elfasta_ready = Channel.empty()
    def elprep_used = callers.contains("elprep")
    if (!ch_elfasta && elprep_used) {
        ELPREP_FASTATOELFASTA(
            ch_fasta_ready
        )
        ch_versions = ch_versions.mix(ELPREP_FASTATOELFASTA.out.versions)
        ch_elfasta_ready = ELPREP_FASTATOELFASTA.out.elfasta
    } else {
        ch_elfasta_ready = ch_elfasta
    }

    // Reference STR table file
    def ch_strtablefile_ready = Channel.empty()
    if (dragstr && !ch_strtablefile) {
        COMPOSESTRTABLEFILE(
            ch_fasta_ready,
            ch_fai_ready,
            ch_dict_ready
        )
        ch_versions  = ch_versions.mix(COMPOSESTRTABLEFILE.out.versions)
        ch_strtablefile_ready = COMPOSESTRTABLEFILE.out.str_table.collect()
    } else if (dragstr) {
        ch_strtablefile_ready = ch_strtablefile
    } else {
        ch_strtablefile_ready = []
    }

    // Reference validation SDF
    def ch_sdf_ready = Channel.empty()
    if (validate && !ch_sdf) {
        RTGTOOLS_FORMAT(
            ch_fasta_ready.map { meta, fasta_file -> [meta, fasta_file, [], []] }
        )
        ch_versions  = ch_versions.mix(RTGTOOLS_FORMAT.out.versions)
        ch_sdf_ready = RTGTOOLS_FORMAT.out.sdf.collect()
    }
    else if (validate && sdf.endsWith(".tar.gz")) {
        UNTAR(
            ch_sdf
        )
        ch_versions = ch_versions.mix(UNTAR.out.versions)

        ch_sdf_ready = UNTAR.out.untar.collect()
    } else if(validate) {
        ch_sdf_ready = ch_sdf
    } else {
        ch_sdf_ready = [[],[]]
    }

    // VEP annotation cache
    def ch_vep_cache_ready = Channel.empty()
    if (!ch_vep_cache && annotate) {
        ENSEMBLVEP_DOWNLOAD(
            Channel.of([[id:"vep_cache"], genome == "hg38" ? "GRCh38" : genome, species, vep_cache_version]).collect()
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

        ch_vep_cache_ready = ENSEMBLVEP_DOWNLOAD.out.cache.collect{ _meta, cache -> cache }
    } else {
        ch_vep_cache_ready = ch_vep_cache
    }

    //
    // Split the input channel into the right channels
    //

    def ch_input = ch_samplesheet
        .multiMap { meta, cram, crai, gvcf, tbi, roi_file, truth_vcf, truth_tbi, truth_bed ->
            // Divide the input files into their corresponding channel
            def new_meta = meta + [
                type: gvcf && cram ? "gvcf_cram" : gvcf ? "gvcf" : "cram" // Define the type of input data
            ]

            def new_meta_validation = meta.subMap(["id", "sample", "family", "duplicate_count"])

            truth_variants: [new_meta_validation, truth_vcf, truth_tbi, truth_bed] // Optional channel containing the truth VCF, its index and the optional BED file
            gvcf:           [new_meta, gvcf, tbi] // Optional channel containing the GVCFs and their optional indices
            cram:           [new_meta, cram, crai]  // Mandatory channel containing the CRAM files and their optional indices
            roi:            [new_meta, roi_file] // Optional channel containing the ROI BED files for WES samples
        }

    //
    // Create the GVCF index if it's missing
    //

    def ch_gvcf_branch = ch_input.gvcf
        .filter { meta, _gvcf, _tbi ->
            // Filter out samples that have no GVCF
            meta.type == "gvcf" || meta.type == "gvcf_cram"
        }
        .branch { meta, gvcf, tbi ->
            no_tbi: !tbi
                return [ meta, gvcf ]
            tbi:    tbi
                return [ meta, gvcf, tbi ]
        }

    TABIX_GVCF(
        ch_gvcf_branch.no_tbi
    )
    ch_versions = ch_versions.mix(TABIX_GVCF.out.versions)

    def ch_gvcfs_ready = ch_gvcf_branch.no_tbi
        .join(TABIX_GVCF.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_gvcf_branch.tbi)
        .combine(callers.intersect(GlobalVariables.gvcfCallers))
        .map { meta, gvcf, tbi, caller ->
            def new_meta = meta + [caller:caller]
            [ new_meta, gvcf, tbi ]
        }

    //
    // Run sample preparation
    //

    def create_bam_files = callers.intersect(GlobalVariables.bamCallers).size() > 0 // Only create BAM files when needed
    CRAM_PREPARE_SAMTOOLS_BEDTOOLS(
        ch_input.cram.filter { meta, _cram, _crai ->
            // Filter out files that already have a called GVCF when only GVCF callers are used
            meta.type == "cram" || (meta.type == "gvcf_cram" && callers - GlobalVariables.gvcfCallers)
        },
        ch_input.roi.filter { meta, _roi_file ->
            // Filter out files that already have a called GVCF when only GVCF callers are used
            meta.type == "cram" || (meta.type == "gvcf_cram" && callers - GlobalVariables.gvcfCallers)
        },
        ch_fasta_ready,
        ch_fai_ready,
        ch_default_roi,
        create_bam_files
    )
    ch_versions = ch_versions.mix(CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.versions)

    //
    // Split the BED files
    //

    def ch_split_cram_bam = Channel.empty()
    if(create_bam_files) {
        ch_split_cram_bam = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams
            .join(CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_bams, failOnDuplicate:true, failOnMismatch:true)
    } else {
        ch_split_cram_bam = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams
    }

    INPUT_SPLIT_BEDTOOLS(
        CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_beds.map { meta, bed ->
            [meta, bed, scatter_count]
        },
        ch_split_cram_bam
    )
    ch_versions = ch_versions.mix(INPUT_SPLIT_BEDTOOLS.out.versions)

    def ch_caller_inputs = INPUT_SPLIT_BEDTOOLS.out.split
        .multiMap { meta, cram, crai, bam=[], bai=[], bed ->
            cram: [meta, cram, crai, bed]
            bam: [meta, bam, bai, bed]
        }

    def ch_calls = Channel.empty()
    if("haplotypecaller" in callers) {
        //
        // Call variants with GATK4 HaplotypeCaller
        //

        CRAM_CALL_GATK4(
            ch_caller_inputs.cram.filter { meta, _cram, _crai, _bed ->
                // Filter out the entries that already have a GVCF
                meta.type == "cram"
            },
            ch_fasta_ready,
            ch_fai_ready,
            ch_dict_ready,
            ch_strtablefile_ready,
            ch_dbsnp_ready,
            ch_dbsnp_tbi_ready,
            dragstr
        )
        ch_gvcfs_ready = ch_gvcfs_ready.mix(CRAM_CALL_GATK4.out.gvcfs)
        ch_versions = ch_versions.mix(CRAM_CALL_GATK4.out.versions)
        ch_reports  = ch_reports.mix(CRAM_CALL_GATK4.out.reports)
    }

    if("elprep" in callers) {
        //
        // Call variants with Elprep
        //

        BAM_CALL_ELPREP(
            ch_caller_inputs.bam.filter { meta, _bam, _bai, _bed ->
                // Filter out the entries that already have a GVCF
                meta.type == "cram"
            },
            ch_elfasta_ready,
            ch_elsites,
            ch_dbsnp_ready,
            ch_dbsnp_tbi_ready
        )
        ch_gvcfs_ready = ch_gvcfs_ready.mix(BAM_CALL_ELPREP.out.gvcfs)
        ch_versions = ch_versions.mix(BAM_CALL_ELPREP.out.versions)
        ch_reports  = ch_reports.mix(BAM_CALL_ELPREP.out.reports)

    }

    if("vardict" in callers) {
        //
        // Call variants with VarDict
        //

        BAM_CALL_VARDICTJAVA(
            ch_caller_inputs.bam,
            ch_fasta_ready,
            ch_fai_ready,
            ch_dbsnp_ready,
            ch_dbsnp_tbi_ready,
            filter
        )
        ch_versions = ch_versions.mix(BAM_CALL_VARDICTJAVA.out.versions)

        ch_calls = ch_calls.mix(BAM_CALL_VARDICTJAVA.out.vcfs)
    }

    // Stop pipeline execution when only calls should happen
    def ch_gvcfs_final = ch_gvcfs_ready.filter { !only_call }

    GVCF_JOINT_GENOTYPE_GATK4(
        ch_gvcfs_final,
        ch_fasta_ready,
        ch_fai_ready,
        ch_dict_ready,
        ch_dbsnp_ready,
        ch_dbsnp_tbi_ready,
        only_merge,
        scatter_count
    )
    ch_versions = ch_versions.mix(GVCF_JOINT_GENOTYPE_GATK4.out.versions)
    ch_calls = ch_calls.mix(GVCF_JOINT_GENOTYPE_GATK4.out.vcfs)

    // Stop pipeline execution when only the merge should happen
    def ch_calls_final = ch_calls.filter { !only_merge }

    def ch_called_variants = ch_calls_final
        .map { meta, vcf, tbi ->
            def new_meta = meta - meta.subMap(["type", "vardict_min_af"])
            [ new_meta, vcf, tbi ]
        }

    BCFTOOLS_STATS(
        ch_called_variants,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
    ch_reports = ch_reports.mix(BCFTOOLS_STATS.out.stats.collect { _meta, report -> report })

    def ch_filtered_variants = Channel.empty()
    if(filter) {
        VCF_FILTER_BCFTOOLS(
            ch_called_variants,
            true
        )
        ch_versions = ch_versions.mix(VCF_FILTER_BCFTOOLS.out.versions)
        ch_filtered_variants = VCF_FILTER_BCFTOOLS.out.vcfs
    } else {
        ch_filtered_variants = ch_called_variants
    }

    def ch_normalized_variants = Channel.empty()
    if(normalize) {
        BCFTOOLS_NORM(
            ch_filtered_variants,
            ch_fasta_ready,
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

        TABIX_NORMALIZE(
            BCFTOOLS_NORM.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_NORMALIZE.out.versions.first())

        ch_normalized_variants = BCFTOOLS_NORM.out.vcf
            .join(TABIX_NORMALIZE.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    } else {
        ch_normalized_variants = ch_filtered_variants
    }

    if(!only_merge && !only_call) {

        //
        // Preprocess the PED channel
        //

        def ch_somalier_input = ch_normalized_variants
            .map { meta, _vcf, _tbi ->
                [ meta, pedFiles.containsKey(meta.family) ? pedFiles[meta.family] : [] ]
            }

        //
        // Run relation tests with somalier
        //

        VCF_EXTRACT_RELATE_SOMALIER(
            ch_normalized_variants,
            ch_fasta_ready,
            ch_fai_ready,
            ch_somalier_sites,
            ch_somalier_input
        )
        ch_versions = ch_versions.mix(VCF_EXTRACT_RELATE_SOMALIER.out.versions)

        //
        // Add PED headers to the VCFs
        //

        def ch_ped_vcfs = Channel.empty()
        if(add_ped){

            VCF_PED_RTGTOOLS(
                ch_normalized_variants,
                VCF_EXTRACT_RELATE_SOMALIER.out.peds
            )
            ch_versions = ch_versions.mix(VCF_PED_RTGTOOLS.out.versions)

            ch_ped_vcfs = VCF_PED_RTGTOOLS.out.ped_vcfs
        } else {
            ch_ped_vcfs = ch_normalized_variants
                .map { meta, vcf, _tbi=[] ->
                    [ meta, vcf ]
                }
        }

        //
        // Annotation of the variants and creation of Gemini-compatible database files
        //

        def ch_annotation_output = Channel.empty()
        if (annotate) {
            VCF_ANNOTATION(
                ch_ped_vcfs,
                ch_fasta_ready,
                ch_vep_cache_ready,
                ch_vep_extra_files,
                ch_vcfanno_config,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                genome,
                species,
                vep_cache_version,
                vep_chunk_size,
                vcfanno
            )
            ch_versions = ch_versions.mix(VCF_ANNOTATION.out.versions)
            ch_reports  = ch_reports.mix(VCF_ANNOTATION.out.reports)

            ch_annotation_output = VCF_ANNOTATION.out.annotated_vcfs
        } else {
            ch_annotation_output = ch_ped_vcfs
        }

        //
        // Tabix the resulting VCF
        //

        TABIX_FINAL(
            ch_annotation_output
        )
        ch_versions = ch_versions.mix(TABIX_FINAL.out.versions.first())

        def ch_final_vcfs = ch_annotation_output
            .join(TABIX_FINAL.out.tbi, failOnDuplicate:true, failOnMismatch:true)

        //
        // Validate the found variants
        //

        if (validate){

            def ch_truths_input = ch_input.truth_variants
                .map { meta, vcf, tbi, bed ->
                    def new_meta = meta - meta.subMap("duplicate_count")
                    [ groupKey(new_meta, meta.duplicate_count), vcf, tbi, bed ]
                }
                .groupTuple()
                .map { meta, vcf, tbi, bed ->
                    // Get only one VCF for samples that were given multiple times
                    def one_vcf = vcf.find { vcf_file -> vcf_file != [] } ?: []
                    def one_tbi = tbi.find { tbi_file -> tbi_file != [] } ?: []
                    def one_bed = bed.find { bed_file -> bed_file != [] } ?: []
                    [ meta, one_vcf, one_tbi, one_bed ]
                }
                .branch { _meta, vcf, tbi, _bed ->
                    no_vcf: !vcf
                    tbi: tbi
                    no_tbi: !tbi
                }

            // Create truth VCF indices if none were given
            TABIX_TRUTH(
                ch_truths_input.no_tbi.map { meta, vcf, _tbi, _bed ->
                    [ meta, vcf ]
                }
            )
            ch_versions = ch_versions.mix(TABIX_TRUTH.out.versions.first())

            ch_truths_input.no_tbi
                .join(TABIX_TRUTH.out.tbi, failOnDuplicate:true, failOnMismatch:true)
                .map { meta, vcf, _empty, bed, tbi ->
                    [ meta, vcf, tbi, bed ]
                }
                .mix(ch_truths_input.tbi)
                .mix(ch_truths_input.no_vcf)
                .combine(callers)
                .map { meta, vcf, tbi, bed, caller ->
                    def new_meta = meta + [caller: caller]
                    [ new_meta, vcf, tbi, bed ]
                }
                .set { ch_truths } // Set needs to be used here due to some Nextflow bug

            def ch_validation_input = ch_final_vcfs
                .map { meta, vcf, tbi ->
                    def new_meta = meta - meta.subMap("family_samples")
                    [ new_meta, vcf, tbi, meta.family_samples.tokenize(",") ]
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
                .join(ch_truths, failOnMismatch:true, failOnDuplicate:true)
                .filter { _meta, _vcf, _tbi, truth_vcf, _truth_tbi, _truth_bed ->
                    // Filter out all samples that have no truth VCF
                    truth_vcf != []
                }
                .multiMap { meta, vcf, tbi, truth_vcf, truth_tbi, truth_bed ->
                    vcfs: [meta, vcf, tbi, truth_vcf, truth_tbi]
                    bed:  [meta, truth_bed]
                }

            CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_beds
                .combine(callers)
                .map { meta, bed, caller ->
                    def new_meta = [
                        id:meta.id,
                        sample:meta.sample,
                        family:meta.family,
                        caller:caller
                    ]
                    [ new_meta, bed ]
                }
                .join(ch_validation_input.bed, failOnMismatch:true, failOnDuplicate:true)
                .map { meta, regions, truth ->
                    [ meta, truth, regions ]
                }
                .set { ch_validation_regions } // Set needs to be used here due to some Nextflow bug

            VCF_VALIDATE_SMALL_VARIANTS(
                ch_validation_input.vcfs,
                ch_validation_regions,
                ch_sdf_ready.collect()
            )
            ch_versions = ch_versions.mix(VCF_VALIDATE_SMALL_VARIANTS.out.versions)
        }

        //
        // Create Gemini-compatible database files
        //

        if(gemini){
            def ch_vcf2db_input = CustomChannelOperators.joinOnKeys(
                    ch_final_vcfs.map { meta, vcf, _tbi -> [ meta, vcf ]},
                    VCF_EXTRACT_RELATE_SOMALIER.out.peds,
                    ['id', 'family', 'family_samples']
                )

            VCF2DB(
                ch_vcf2db_input
            )
            ch_versions = ch_versions.mix(VCF2DB.out.versions.first())

        }

        //
        // Run UPDio analysis
        //

        if(updio) {
            VCF_UPD_UPDIO(
                ch_final_vcfs,
                VCF_EXTRACT_RELATE_SOMALIER.out.peds,
                ch_updio_common_cnvs
            )
            ch_versions = ch_versions.mix(VCF_UPD_UPDIO.out.versions.first())
        }

        //
        // Run automap analysis
        //

        if(automap) {
            VCF_ROH_AUTOMAP(
                ch_final_vcfs,
                ch_automap_repeats,
                ch_automap_panel,
                genome
            )
            ch_versions = ch_versions.mix(VCF_ROH_AUTOMAP.out.versions.first())
        }
    }

    //
    // Collate and save software versions
    //
    def ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // Perform multiQC on all QC data
    //

    def ch_multiqc_config                     = Channel.fromPath(
                                                "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config              = multiqc_config ?
                                                Channel.fromPath(multiqc_config, checkIfExists: true) :
                                                Channel.empty()
    def ch_multiqc_logo                       = multiqc_logo ?
                                                Channel.fromPath(multiqc_logo, checkIfExists: true) :
                                                Channel.empty()

    def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = multiqc_methods_description ?
                                                file(multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files                          = ch_multiqc_files.mix(
                                                    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                                                    ch_collated_versions,
                                                    ch_methods_description.collectFile(
                                                        name: 'methods_description_mqc.yaml',
                                                        sort: false
                                                    ),
                                                    ch_reports
                                                )


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
