/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_cmgg_smallvariants_pipeline'

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
include { CRAM_REPEAT_EXPANSIONHUNTER       } from '../subworkflows/local/cram_repeat_expansionhunter/main'
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
include { BCFTOOLS_QUERY as BCFTOOLS_GETSAMPLES                      } from '../modules/nf-core/bcftools/query/main'
include { TABIX_TABIX as TABIX_DBSNP                                 } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GVCF                                  } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VCF                                   } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TRUTH                                 } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_FAMILY                    } from '../modules/nf-core/bcftools/stats/main'
include { VCF2DB                                                     } from '../modules/nf-core/vcf2db/main'
include { MULTIQC                                                    } from '../modules/nf-core/multiqc/main'
include { MSISENSORPRO_PRO                                           } from '../modules/nf-core/msisensorpro/pro/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The main workflow
workflow SMALLVARIANTS {

    take:
    // Input channels
    ch_samplesheet              // dataflow queue: The input channel
    ped_files                   // dataflow value: A channel with a map of family keys and PED files as values

    // File inputs
    fasta                       // path: path to the reference fasta
    fai                         // path: path to the index of the reference fasta
    dict                        // path: path to the sequence dictionary file
    elfasta                     // path: path to the elfasta reference file
    strtablefile                // path: path to the strtable file
    sdf                         // path: path to the SDF directory
    dbsnp                       // path: path to the DBSNP VCF file
    dbsnp_tbi                   // path: path to the index of the DBSNP VCF file
    vep_cache                   // path: path to the VEP cache
    dbnsfp                      // path: path to the DBNSFP file
    dbnsfp_tbi                  // path: path to the index of the DBNSFP file
    spliceai_indel              // path: path to the SpliceAI indels file
    spliceai_indel_tbi          // path: path to the index of the SpliceAI indels file
    spliceai_snv                // path: path to the SpliceAI SNV file
    spliceai_snv_tbi            // path: path to the index of the SpliceAI SNV file
    mastermind                  // path: path to the Mastermind file
    mastermind_tbi              // path: path to the index of the Mastermind file
    eog                         // path: path to the EOG file
    eog_tbi                     // path: path to the index of the EOG file
    alphamissense               // path: path to the Alphamissense file
    alphamissense_tbi           // path: path to the index of the Alphamissense file
    maxentscan                  // path: path to the MaxEntScan directory
    vcfanno_resources           // string: semicolon-separated paths/globs to the VCFanno resources
    vcfanno_config              // path: path to VCFanno config TOML
    multiqc_config              // path: path to the multiqc config file
    multiqc_logo                // path: path to the multiqc logo
    multiqc_methods_description // path: path to the multiqc methods description
    roi                         // path: path to the default ROI BED file
    somalier_sites              // path: path to the Somalier sites file
    vcfanno_lua                 // path: path to the VCFanno Lua script
    updio_common_cnvs           // path: path to the file containing common UPDio CNVs
    automap_repeats             // path: path to the Automap repeats file
    automap_panel               // path: path to the Automap panel file
    outdir                      // path: path to the output directory
    elsites                     // path: path to the elsites file for elprep
    msi_baseline                // path: path to the msi_baseline file
    updio_regions               // path: path to the BED file with regions to be used by UPDio
    expansionhunter_catalogue     // path: path to the ExpansionHunter variant catalog

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
    disable_updio               // boolean: disable running UPDio on the final VCFs
    disable_automap             // boolean: disable running Automap on the final VCFs
    vep_dbnsfp                  // boolean: use the DBNSFP VEP plugin
    vep_spliceai                // boolean: use the SpliceAI VEP plugin
    vep_mastermind              // boolean: use the Mastermind VEP plugin
    vep_eog                     // boolean: use the EOG VEP plugin
    vep_alphamissense           // boolean: use the AlphaMissense VEP plugin
    vep_maxentscan              // boolean: use the MaxEntScan VEP plugin

    // Value inputs
    genome                      // string:  the genome used by the pipeline run
    species                     // string:  the species used by the pipeline run
    vep_cache_version           // integer: the vep cache version to be used
    vep_chunk_size              // integer: the chunk size to split each VCF file into for VEP
    scatter_count               // integer: the amount of scattering performed on each file
    callers                     // list:    the callers to use


    main:
    def ch_reports       = channel.empty()
    def ch_multiqc_files = channel.empty()

    def List<String> gvcf_callers = ["haplotypecaller", "elprep"]
    def List<String> bam_callers = ["elprep", "vardict"]
    def List<String> snv_callers = ["haplotypecaller", "vardict", "elprep"]

    //
    // Importing and convert the input files passed through the parameters to channels
    //

    def ch_fasta_ready        = channel.value([[id:"reference"], fasta])
    def ch_fai                = fai                 ? channel.value([[id:"reference"], fai]) : null
    def ch_dict               = dict                ? channel.value([[id:"reference"], dict]) : null
    def ch_elfasta            = elfasta             ? channel.value([[id:"reference"], elfasta]) : null
    def ch_strtablefile       = strtablefile        ? channel.value([[id:"reference"], strtablefile]) : null
    def ch_sdf                = sdf                 ? channel.value([[id:'reference'], sdf]) : null

    def ch_default_roi        = roi                 ? channel.value(roi) : []

    def ch_dbsnp_ready        = dbsnp               ? channel.value([[id:"dbsnp"], dbsnp]) : [[],[]]
    def ch_dbsnp_tbi          = dbsnp_tbi           ? channel.value([[id:"dbsnp"], dbsnp_tbi]) : [[],[]]

    def ch_somalier_sites     = somalier_sites      ? channel.value([[id:"somalier_sites"], somalier_sites]) : [[],[]]

    def ch_vep_cache          = vep_cache           ? channel.value([[id:'cache'], file(vep_cache)]) : []

    def ch_vcfanno_config     = vcfanno_config      ? channel.value(vcfanno_config) : []
    def ch_vcfanno_lua        = vcfanno_lua         ? channel.value(vcfanno_lua) : []
    def ch_vcfanno_resources  = vcfanno_resources   ? channel.value(vcfanno_resources.split(";").collect{ res -> files(res, checkIfExists:true) }.flatten()) : []

    def ch_updio_common_cnvs  = updio_common_cnvs   ? channel.value([[id:'updio_cnv'], updio_common_cnvs]) : [[],[]]

    def ch_automap_repeats    = automap_repeats     ? channel.value([[id:"repeats"], automap_repeats]) : []
    def ch_automap_panel      = automap_panel       ? channel.value([[id:"automap_panel"], automap_panel]) : [[],[]]

    def ch_elsites            = elsites             ? channel.fromPath(elsites).map{ elsites_file -> [[id:'elsites'], elsites_file] }.collect() : [[],[]]

    def ch_msi_baseline       = msi_baseline        ? channel.value([[id:"msi_baseline"], msi_baseline]) : [[],[]]

    def ch_updio_regions      = updio_regions       ? channel.value(updio_regions) : []

    def ch_expansionhunter_catalogue = expansionhunter_catalogue ? channel.value([[id:"expansionhunter_catalogue"], expansionhunter_catalogue]) : channel.empty()

    //
    // Check for the presence of EnsemblVEP plugins that use extra files
    //

    def ch_vep_extra_files = []

    if(annotate){
        // Check if all dbnsfp files are given
        if (dbnsfp && dbnsfp_tbi && vep_dbnsfp) {
            ch_vep_extra_files.add(dbnsfp)
            ch_vep_extra_files.add(dbnsfp_tbi)
        }
        else if (vep_dbnsfp) {
            error("Please specify '--vep_dbsnfp true', '--dbnsfp PATH/TO/DBNSFP/FILE' and '--dbnspf_tbi PATH/TO/DBNSFP/INDEX/FILE' to use the dbnsfp VEP plugin.")
        }

        // Check if all spliceai files are given
        if (spliceai_snv && spliceai_snv_tbi && spliceai_indel && spliceai_indel_tbi && vep_spliceai) {
            ch_vep_extra_files.add(spliceai_snv)
            ch_vep_extra_files.add(spliceai_snv_tbi)
            ch_vep_extra_files.add(spliceai_indel)
            ch_vep_extra_files.add(spliceai_indel_tbi)
        }
        else if (vep_spliceai) {
            error("Please specify '--vep_spliceai true', '--spliceai_snv PATH/TO/SPLICEAI/SNV/FILE', '--spliceai_snv_tbi PATH/TO/SPLICEAI/SNV/INDEX/FILE', '--spliceai_indel PATH/TO/SPLICEAI/INDEL/FILE' and '--spliceai_indel_tbi PATH/TO/SPLICEAI/INDEL/INDEX/FILE' to use the SpliceAI VEP plugin.")
        }

        // Check if all mastermind files are given
        if (mastermind && mastermind_tbi && vep_mastermind) {
            ch_vep_extra_files.add(mastermind)
            ch_vep_extra_files.add(mastermind_tbi)
        }
        else if (vep_mastermind) {
            error("Please specify '--vep_mastermind true', '--mastermind PATH/TO/MASTERMIND/FILE' and '--mastermind_tbi PATH/TO/MASTERMIND/INDEX/FILE' to use the mastermind VEP plugin.")
        }

        // Check if all EOG files are given
        if (eog && eog_tbi && vep_eog) {
            ch_vep_extra_files.add(eog)
            ch_vep_extra_files.add(eog_tbi)
        }
        else if (vep_eog) {
            error("Please specify '--vep_eog true', '--eog PATH/TO/EOG/FILE' and '--eog_tbi PATH/TO/EOG/INDEX/FILE' to use the EOG custom VEP plugin.")
        }

        // Check if all AlphaMissense files are given
        if (alphamissense && alphamissense_tbi && vep_alphamissense) {
            ch_vep_extra_files.add(alphamissense)
            ch_vep_extra_files.add(alphamissense_tbi)
        }
        else if (vep_alphamissense) {
            error("Please specify '--vep_alphamissense true', '--alphamissense PATH/TO/ALPHAMISSENSE/FILE' and '--alphamissense_tbi PATH/TO/ALPHAMISSENSE/INDEX/FILE' to use the AlphaMissense VEP plugin.")
        }

        // Check if all maxentscan files are given
        if (maxentscan && vep_maxentscan) {
            ch_vep_extra_files.add(maxentscan)
        }
        else if (vep_maxentscan) {
            error("Please specify '--vep_maxentscan true' and '--maxentscan PATH/TO/MAXENTSCAN/DIRECTORY' to use the MaxEntScan VEP plugin.")
        }
    }

    //
    // Create the optional input files if they are not supplied
    //

    // DBSNP index
    def ch_dbsnp_tbi_ready = channel.empty()
    if (dbsnp && !dbsnp_tbi) {
        TABIX_DBSNP(
            ch_dbsnp_ready.map { meta, dbsnp_file -> [meta, dbsnp_file, [], []] }
        )
        ch_dbsnp_tbi_ready = TABIX_DBSNP.out.index.collect()
    } else {
        ch_dbsnp_tbi_ready = ch_dbsnp_tbi
    }

    // Reference fasta index
    def ch_fai_ready = channel.empty()
    if (!fai) {
        FAIDX(
            ch_fasta_ready,
            [[],[]]
        )
        ch_fai_ready = FAIDX.out.fai
            .collect()
    } else {
        ch_fai_ready = ch_fai
    }

    // Reference sequence dictionary
    def ch_dict_ready = channel.empty()
    if (!dict) {
        CREATESEQUENCEDICTIONARY(
            ch_fasta_ready
        )
        ch_dict_ready = CREATESEQUENCEDICTIONARY.out.dict.collect()
    } else {
        ch_dict_ready = ch_dict
    }

    def ch_elfasta_ready = channel.empty()
    def elprep_used = callers.contains("elprep")
    if (!elfasta && elprep_used) {
        ELPREP_FASTATOELFASTA(
            ch_fasta_ready
        )
        ch_elfasta_ready = ELPREP_FASTATOELFASTA.out.elfasta
    } else {
        ch_elfasta_ready = ch_elfasta
    }

    // Reference STR table file
    def ch_strtablefile_ready = channel.empty()
    if (dragstr && !strtablefile) {
        COMPOSESTRTABLEFILE(
            ch_fasta_ready,
            ch_fai_ready,
            ch_dict_ready
        )
        ch_strtablefile_ready = COMPOSESTRTABLEFILE.out.str_table.collect()
    } else if (dragstr) {
        ch_strtablefile_ready = ch_strtablefile
    } else {
        ch_strtablefile_ready = []
    }

    // Reference validation SDF
    def ch_sdf_ready = channel.empty()
    if (validate && !sdf) {
        RTGTOOLS_FORMAT(
            ch_fasta_ready.map { meta, fasta_file -> [meta, fasta_file, [], []] }
        )
        ch_sdf_ready = RTGTOOLS_FORMAT.out.sdf.collect()
    }
    else if (validate && sdf.toUriString().endsWith(".tar.gz")) {
        UNTAR(
            ch_sdf
        )
        ch_sdf_ready = UNTAR.out.untar.collect()
    } else if(validate) {
        ch_sdf_ready = ch_sdf
    } else {
        ch_sdf_ready = [[],[]]
    }

    // VEP annotation cache
    def ch_vep_cache_ready = channel.empty()
    if (!vep_cache && annotate) {
        ENSEMBLVEP_DOWNLOAD(
            channel.of([[id:"vep_cache"], genome == "hg38" ? "GRCh38" : genome, species, vep_cache_version]).collect(),
            false
        )
        ch_vep_cache_ready = ENSEMBLVEP_DOWNLOAD.out.cache.collect()
    } else {
        ch_vep_cache_ready = ch_vep_cache
    }

    //
    // Split the input channel into the right channels
    //

    def usedGvcfCallers = callers.intersect(gvcf_callers)

    def ch_input = ch_samplesheet
        .multiMap { meta, cram, crai, gvcf, gtbi, vcf, tbi, roi_file, truth_vcf, truth_tbi, truth_bed ->
            // Error checks that were not possible using nf-schema
            if (gvcf && usedGvcfCallers.size() >= 2) {
                error("GVCF input is not supported for runs that use more than one caller that produces a GVCF output")
            }
            if (gvcf && validate) {
                error("Validation is not supported for GVCF inputs, use CRAM files instead when using `--validate`.")
            }
            if (vcf && validate) {
                error("Validation is not supported for VCF inputs, use CRAM files instead when using `--validate`.")
            }

            // Don't analyse CRAM and GVCF when VCF is given
            if (vcf && (cram || crai || gvcf || gtbi)) {
                log.warn("Found a VCF file for family ${meta.family}, skipping all analysis on CRAM and GVCF files")
                cram = []
                crai = []
                gvcf = []
                gtbi = []
            }

            // Divide the input files into their corresponding channel
            def new_meta = meta + [
                type: vcf ? "vcf" : gvcf && cram ? "gvcf_cram" : gvcf ? "gvcf" : "cram" // Define the type of input data
            ]

            def new_meta_validation = meta.subMap(["id", "sample", "family", "duplicate_count"])
            def new_meta_vcf = meta + [id:meta.family, caller:"unknown_caller"] - meta.subMap(["sample", "vardict_min_af"])

            def new_meta_gvcf = meta
            if (usedGvcfCallers.size() == 1) {
                new_meta_gvcf = meta + [caller:usedGvcfCallers[0]]
            }

            truth_variants: [new_meta_validation, truth_vcf, truth_tbi, truth_bed] // Optional channel containing the truth VCF, its index and the optional BED file
            gvcf:           [new_meta, gvcf, gtbi] // Optional channel containing the GVCFs and their optional indices
            cram:           [new_meta, cram, crai]  // Mandatory channel containing the CRAM files and their optional indices
            roi:            [new_meta, roi_file] // Optional channel containing the ROI BED files for WES samples
            vcf:            [new_meta_vcf, vcf, tbi] // Optional channel containing the VCFs and their optional indices
        }

    //
    // Handle the vcf branch
    //

    def ch_vcf_branch = ch_input.vcf
        .filter { _meta, vcf, _tbi -> vcf } // Filter out samples that have no VCF
        .branch { meta, vcf, tbi ->
            no_tbi: !tbi
                return [ meta, vcf ]
            tbi:    tbi
                return [ meta, vcf, tbi ]
        }

    TABIX_VCF(
        ch_vcf_branch.no_tbi.map { meta, vcf -> [meta, vcf, [], []] }
    )

    def ch_indexed_vcfs = ch_vcf_branch.no_tbi
        .join(TABIX_VCF.out.index, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_vcf_branch.tbi)

    BCFTOOLS_GETSAMPLES(
        ch_indexed_vcfs,
        [],
        [],
        []
    )

    def ch_vcfs_ready = BCFTOOLS_GETSAMPLES.out.output
        .join(ch_indexed_vcfs, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, samples, vcf, tbi ->
            def new_meta = meta + [family_samples: samples.text.readLines().join(",")]
            [ new_meta, vcf, tbi ]
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
        ch_gvcf_branch.no_tbi.map { meta, gvcf -> [meta, gvcf, [], []] }
    )

    def ch_gvcfs_ready = ch_gvcf_branch.no_tbi
        .join(TABIX_GVCF.out.index, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_gvcf_branch.tbi)
        .combine(callers.intersect(gvcf_callers))
        .map { meta, gvcf, tbi, caller ->
            def new_meta = meta + [caller:caller]
            [ new_meta, gvcf, tbi ]
        }

    //
    // Run sample preparation
    //

    def create_bam_files = callers.intersect(bam_callers).size() > 0 // Only create BAM files when needed
    CRAM_PREPARE_SAMTOOLS_BEDTOOLS(
        ch_input.cram.filter { meta, _cram, _crai ->
            // Filter out files that already have a called GVCF when only GVCF callers are used
            meta.type == "cram" || (meta.type == "gvcf_cram" && callers - gvcf_callers)
        },
        ch_input.roi.filter { meta, _roi_file ->
            // Filter out files that already have a called GVCF when only GVCF callers are used
            meta.type == "cram" || (meta.type == "gvcf_cram" && callers - gvcf_callers)
        },
        ch_fasta_ready,
        ch_fai_ready,
        ch_default_roi,
        create_bam_files
    )
    ch_reports  = ch_reports.mix(CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.reports)
    def ch_single_beds = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_beds
    def ch_perbase_beds = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.perbase_beds
    def ch_merged_crams = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.merged_crams
    def ch_mosdepth_reports = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.mosdepth_reports

    //
    // Call repeat expansions
    //

    def ch_expansionhunter_vcfs = channel.empty()
    if ("expansionhunter" in callers) {
        CRAM_REPEAT_EXPANSIONHUNTER(
            CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams,
            ch_fasta_ready,
            ch_fai_ready,
            ch_expansionhunter_catalogue
        )
        ch_expansionhunter_vcfs = CRAM_REPEAT_EXPANSIONHUNTER.out
    }

    //
    // Split the BED files
    //

    def ch_split_cram_bam = channel.empty()
    if(create_bam_files) {
        ch_split_cram_bam = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams
            .join(CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_bams, failOnDuplicate:true, failOnMismatch:true)
    } else {
        ch_split_cram_bam = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams
    }

    INPUT_SPLIT_BEDTOOLS(
        ch_single_beds.map { meta, bed ->
            [meta, bed, scatter_count]
        },
        ch_split_cram_bam
    )

    def ch_caller_inputs = INPUT_SPLIT_BEDTOOLS.out.split
        .multiMap { meta, cram, crai, bam=[], bai=[], bed ->
            cram: [meta, cram, crai, bed]
            bam: [meta, bam, bai, bed]
        }

    //
    // Check for MSI
    //

    def msi_warned = false
    def ch_msi_samples = CRAM_PREPARE_SAMTOOLS_BEDTOOLS.out.ready_crams
        .filter { meta, _cram, _crai ->
            if(!msi_baseline && meta.msi) {
                if(!msi_warned) {
                    log.warn("MSI samples were found, but no MSI baseline file was provided. Please provide a baseline file using the '--msi_baseline' parameter. Skipping MSI analysis...")
                }
                msi_warned = true
                return false
            }
            return meta.msi
        }

    MSISENSORPRO_PRO(
        ch_msi_samples,
        ch_msi_baseline,
        ch_fasta_ready,
        ch_fai_ready
    )
    ch_msisensor_output = MSISENSORPRO_PRO.out.all_msi.mix(
        MSISENSORPRO_PRO.out.summary_msi,
        MSISENSORPRO_PRO.out.dis_msi,
        MSISENSORPRO_PRO.out.unstable_msi
    )
    ch_reports  = ch_reports.mix(MSISENSORPRO_PRO.out.all_msi.map { _meta, file -> file})
    ch_reports  = ch_reports.mix(MSISENSORPRO_PRO.out.summary_msi.map { _meta, file -> file})

    def ch_calls = ch_vcfs_ready
    def ch_gvcf_reports = channel.empty()
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
        ch_reports  = ch_reports.mix(CRAM_CALL_GATK4.out.reports.map { _meta, report -> report })
        ch_gvcf_reports = ch_gvcf_reports.mix(CRAM_CALL_GATK4.out.reports)
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
        ch_reports  = ch_reports.mix(BAM_CALL_ELPREP.out.reports.map { _meta, report -> report })
        ch_gvcf_reports = ch_gvcf_reports.mix(BAM_CALL_ELPREP.out.reports)
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
            ch_dbsnp_tbi_ready
        )
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
    ch_calls = ch_calls.mix(GVCF_JOINT_GENOTYPE_GATK4.out.vcfs)
    def ch_joint_beds = GVCF_JOINT_GENOTYPE_GATK4.out.beds
    def ch_final_genomicsdb = GVCF_JOINT_GENOTYPE_GATK4.out.genomicsdb

    def ch_final_vcfs       = channel.empty()
    def ch_final_dbs        = channel.empty()
    def ch_final_peds       = channel.empty()
    def ch_final_reports    = channel.empty()
    def ch_final_automap    = channel.empty()
    def ch_final_updio      = channel.empty()
    def ch_final_validation = channel.empty()

    if (!only_call && !only_merge) {
        def ch_called_variants = ch_calls
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
        ch_final_reports = BCFTOOLS_STATS.out.stats
        ch_reports = ch_reports.mix(ch_final_reports.collect { _meta, report -> report })

        def ch_filtered_variants = channel.empty()
        if(filter) {
            VCF_FILTER_BCFTOOLS(
                ch_called_variants
            )
            ch_filtered_variants = VCF_FILTER_BCFTOOLS.out.vcfs
        } else {
            ch_filtered_variants = ch_called_variants
        }

        def ch_normalized_variants = channel.empty()
        if(normalize) {
            BCFTOOLS_NORM(
                ch_filtered_variants,
                ch_fasta_ready,
            )
            ch_normalized_variants = BCFTOOLS_NORM.out.vcf
                .join(BCFTOOLS_NORM.out.index, failOnDuplicate:true, failOnMismatch:true)
        } else {
            ch_normalized_variants = ch_filtered_variants
        }

        //
        // Preprocess the PED channel
        //

        def ch_somalier_input = ch_normalized_variants
            .combine(ped_files)
            .map { meta, _vcf, _tbi, ped_map ->
                [ meta, ped_map.get(meta.family, []) ]
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
        ch_final_peds = VCF_EXTRACT_RELATE_SOMALIER.out.peds
        ch_final_reports = ch_final_reports.mix(VCF_EXTRACT_RELATE_SOMALIER.out.html)
        ch_reports = ch_reports.mix(VCF_EXTRACT_RELATE_SOMALIER.out.pairs_tsv.map { _meta, report -> report })
        ch_reports = ch_reports.mix(VCF_EXTRACT_RELATE_SOMALIER.out.samples_tsv.map { _meta, report -> report })

        //
        // Add PED headers to the VCFs
        //

        def ch_ped_vcfs = channel.empty()
        if(add_ped){

            VCF_PED_RTGTOOLS(
                ch_normalized_variants,
                ch_final_peds
            )
            ch_ped_vcfs = VCF_PED_RTGTOOLS.out.ped_vcfs
        } else {
            ch_ped_vcfs = ch_normalized_variants
        }

        //
        // Annotation of the variants and creation of Gemini-compatible database files
        //

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
            ch_reports  = ch_reports.mix(VCF_ANNOTATION.out.reports)

            ch_final_vcfs = VCF_ANNOTATION.out.annotated_vcfs
        } else {
            ch_final_vcfs = ch_ped_vcfs
        }

        //
        // Validate the found variants
        //

        if (validate && callers.intersect(snv_callers)){
            def callers_to_validate = callers.intersect(snv_callers)
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
                    [ meta, vcf, [], [] ]
                }
            )

            ch_truths_input.no_tbi
                .join(TABIX_TRUTH.out.index, failOnDuplicate:true, failOnMismatch:true)
                .map { meta, vcf, _empty, bed, tbi ->
                    [ meta, vcf, tbi, bed ]
                }
                .mix(ch_truths_input.tbi)
                .mix(ch_truths_input.no_vcf)
                .combine(callers_to_validate)
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

            ch_single_beds
                .combine(callers_to_validate)
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

            ch_final_validation = VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_true_positive_vcf.mix(
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_true_positive_vcf_tbi,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_false_negative_vcf,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_false_negative_vcf_tbi,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_false_positive_vcf,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_false_positive_vcf_tbi,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_true_positive_baseline_vcf,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_true_positive_baseline_vcf_tbi,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_summary,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_phasing,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_snp_roc,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_non_snp_roc,
                VCF_VALIDATE_SMALL_VARIANTS.out.vcfeval_weighted_roc,
                VCF_VALIDATE_SMALL_VARIANTS.out.rtgtools_snp_png_rocplot,
                VCF_VALIDATE_SMALL_VARIANTS.out.rtgtools_non_snp_png_rocplot,
                VCF_VALIDATE_SMALL_VARIANTS.out.rtgtools_weighted_png_rocplot,
                VCF_VALIDATE_SMALL_VARIANTS.out.rtgtools_snp_svg_rocplot,
                VCF_VALIDATE_SMALL_VARIANTS.out.rtgtools_non_snp_svg_rocplot,
                VCF_VALIDATE_SMALL_VARIANTS.out.rtgtools_weighted_svg_rocplot
            )
        }

        //
        // Create Gemini-compatible database files
        //

        if(gemini){
            def ch_vcf2db_input = ch_final_vcfs.map { meta, vcf, _tbi -> [ meta, vcf ]}
                .join(ch_final_peds, failOnMismatch:true, failOnDuplicate:true)

            VCF2DB(
                ch_vcf2db_input
            )
            ch_final_dbs = VCF2DB.out.db
        }

        //
        // Run UPDio analysis
        //

        if(!disable_updio) {
            VCF_UPD_UPDIO(
                ch_final_vcfs,
                ch_final_peds,
                ch_updio_common_cnvs,
                ch_updio_regions
            )
            ch_final_updio = VCF_UPD_UPDIO.out.updio
        }

        //
        // Run automap analysis
        //

        if(!disable_automap) {
            VCF_ROH_AUTOMAP(
                ch_final_vcfs,
                ch_automap_repeats,
                ch_automap_panel,
                genome
            )
            ch_final_automap = VCF_ROH_AUTOMAP.out.automap
        }
    }

    //
    // Collate and save software versions
    //
    def versions = channel.topic("versions")
        .distinct()
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_cmgg_smallvariants_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    //
    // Perform multiQC on all QC data
    //

    def ch_multiqc_custom_config              = multiqc_config ?
                                                channel.fromPath(multiqc_config, checkIfExists: true) :
                                                channel.value([])
    def ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
                                                .combine(ch_multiqc_custom_config)
                                                .collect()
                                                .map { configs -> [configs] }
    def ch_multiqc_logo                       = multiqc_logo ?
                                                channel.fromPath(multiqc_logo, checkIfExists: true) :
                                                channel.value([])

    def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = multiqc_methods_description ?
                                                file(multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files                          = ch_multiqc_files.mix(
                                                    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                                                    versions,
                                                    ch_methods_description.collectFile(
                                                        name: 'methods_description_mqc.yaml',
                                                        sort: false
                                                    ),
                                                    ch_reports
                                                )

    def ch_multiqc_input = ch_multiqc_files.collect().map { files -> [files] }
        .combine(ch_multiqc_config)
        .combine(ch_multiqc_logo.toList())
        .map { files, configs, logo ->
            [ [id: 'multiqc'], files, configs, logo, [], [] ]
        }

    MULTIQC(ch_multiqc_input)

    emit:
    merged_crams        = ch_merged_crams               // channel: [ val(meta), path(cram), path(crai) ]
    mosdepth_reports    = ch_mosdepth_reports           // channel: [ val(meta), path(mosdepth_report) ]
    gvcfs               = ch_gvcfs_ready                // channel: [ val(meta), path(gvcf), path(tbi) ]
    msi                 = ch_msisensor_output           // channel: [ val(meta), path(file) ]
    genomicsdb          = ch_final_genomicsdb           // channel: [ val(meta), path(genomicsdb) ]
    vcfs                = ch_final_vcfs                 // channel: [ val(meta), path(vcf), path(tbi) ]
    repeat_vcfs         = ch_expansionhunter_vcfs       // channel: [ val(meta), path(vcf), path(tbi) ]
    gemini              = ch_final_dbs                  // channel: [ val(meta), path(db) ]
    peds                = ch_final_peds                 // channel: [ val(meta), path(ped) ]
    single_beds         = ch_single_beds                // channel: [ val(meta), path(bed) ]
    perbase_beds        = ch_perbase_beds               // channel: [ val(meta), path(bed), path(csi) ]
    joint_beds          = ch_joint_beds                 // channel: [ val(meta), path(bed) ]
    final_reports       = ch_final_reports              // channel: [ val(meta), path(report) ]
    gvcf_reports        = ch_gvcf_reports               // channel: [ val(meta), path(report) ]
    automap             = ch_final_automap              // channel: [ val(meta), path(automap) ]
    updio               = ch_final_updio                // channel: [ val(meta), path(updio) ]
    validation          = ch_final_validation           // channel: [ val(meta), path(file) ]
    multiqc_report      = MULTIQC.out.report            // channel: /path/to/multiqc_report.html
    multiqc_data        = MULTIQC.out.data              // channel: /path/to/multiqc_data
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
