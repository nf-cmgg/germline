//
// SAMPLE_PREPARATION
//

include { MERGE_BEDS as MERGE_ROI_PARAMS    } from '../../../modules/local/merge_beds'
include { MERGE_BEDS as MERGE_ROI_SAMPLE    } from '../../../modules/local/merge_beds'
include { FILTER_BEDS                       } from '../../../modules/local/filter_beds/main'

include { SAMTOOLS_MERGE                    } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT                  } from '../../../modules/nf-core/samtools/convert/main'
include { TABIX_TABIX                       } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP as UNZIP_ROI          } from '../../../modules/nf-core/tabix/bgzip/main'
include { BEDTOOLS_INTERSECT                } from '../../../modules/nf-core/bedtools/intersect/main'
include { MOSDEPTH                          } from '../../../modules/nf-core/mosdepth/main'

workflow CRAM_PREPARE_SAMTOOLS_BEDTOOLS {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their optional indices
        ch_roi               // channel: [optional]  [ val(meta), path(roi) ] => ROI bed files for WES analysis
        ch_fasta             // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_default_roi       // channel: [optional]  [ path(roi) ] => bed containing regions of interest to be used as default
        output_bam           // boolean: Also output BAM files

    main:

    def ch_versions  = Channel.empty()
    def ch_reports   = Channel.empty()

    //
    // Merge the CRAM files if there are multiple per sample
    //

    def ch_cram_branch = ch_crams
        .map { meta, cram, crai ->
            [ groupKey(meta, meta.duplicate_count), cram, crai]
        }
        .groupTuple()
        .branch { meta, cram, crai ->
            multiple: cram.size() > 1
                return [meta.target, cram]
            single:   cram.size() == 1
                return [meta.target, cram[0], crai[0]]
        }

    SAMTOOLS_MERGE(
        ch_cram_branch.multiple,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    //
    // Index the CRAM files which have no index
    //

    def ch_merged_crams = SAMTOOLS_MERGE.out.cram
        .mix(ch_cram_branch.single)
        .branch { meta, cram, crai=[] ->
            not_indexed: crai == []
                return [ meta, cram ]
            indexed: crai != []
                return [ meta, cram, crai ]
        }

    SAMTOOLS_INDEX(
        ch_merged_crams.not_indexed
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def ch_ready_crams = ch_merged_crams.not_indexed
        .join(SAMTOOLS_INDEX.out.crai, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_merged_crams.indexed)

    //
    // Optionally convert the CRAM files to BAM
    //

    def ch_ready_bams = Channel.empty()
    if(output_bam) {
        SAMTOOLS_CONVERT(
            ch_ready_crams,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

        ch_ready_bams = SAMTOOLS_CONVERT.out.bam.join(SAMTOOLS_CONVERT.out.bai, failOnDuplicate:true, failOnMismatch:true)
    }

    //
    // Preprocess the ROI BED files => sort and merge overlapping regions
    //

    def ch_roi_branch = ch_roi
        .map { meta, roi ->
            [ groupKey(meta, meta.duplicate_count), roi ]
        }
        .groupTuple()
        .branch { meta, roi ->
            // Determine whether there is an ROI file given to the current sample
            // It's possible that a sample is given multiple times in the samplesheet, in which
            // case they have been merged earlier. This code checks if at least one entry of the same
            // sample contains an ROI file
            def output_roi = roi.findAll { entry -> entry != [] }
            found:      output_roi.size() > 0
                return [ meta.target, output_roi ]
            missing:    output_roi.size() == 0
                return [ meta.target, [] ]
        }

    MERGE_ROI_SAMPLE(
        ch_roi_branch.found,
        ch_fai
    )
    ch_versions = ch_versions.mix(MERGE_ROI_SAMPLE.out.versions.first())

    // Add the default ROI file to all samples without an ROI file
    // if an ROI BED file has been given through the --roi parameter
    def ch_missing_rois = Channel.empty()
    if (ch_default_roi) {
        MERGE_ROI_PARAMS(
            ch_default_roi.map { bed ->
                [[id:"default_roi"], bed]
            },
            ch_fai
        )
        ch_versions = ch_versions.mix(MERGE_ROI_PARAMS.out.versions)

        ch_missing_rois = ch_roi_branch.missing
            .map { meta, bed ->
                [ groupKey(meta, meta.duplicate_count), bed ]
            }
            .groupTuple()
            .combine(MERGE_ROI_PARAMS.out.bed.map { _meta, bed -> bed })
            .map { meta, _missing, default_roi ->
                [ meta.target, default_roi ]
            }
    } else {
        ch_missing_rois = ch_roi_branch.missing
    }

    def ch_ready_rois = ch_missing_rois.mix(MERGE_ROI_SAMPLE.out.bed)

    //
    // Create callable regions
    //

    def ch_mosdepth_input = ch_ready_crams
        .join(ch_ready_rois, failOnDuplicate:true, failOnMismatch:true)

    MOSDEPTH(
        ch_mosdepth_input,
        ch_fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    def ch_beds_to_filter = ch_ready_rois
        .join(MOSDEPTH.out.quantized_bed, failOnDuplicate:true, failOnMismatch:true)

    // Filter out the regions with no coverage
    FILTER_BEDS(
        ch_beds_to_filter.map { meta, _roi, callable -> [ meta, callable ]}
    )
    ch_versions = ch_versions.mix(FILTER_BEDS.out.versions)

    def ch_beds_to_intersect = FILTER_BEDS.out.bed
        .join(ch_beds_to_filter, failOnDuplicate:true, failOnMismatch:true)
        .branch { meta, filtered_callable, roi, _callable ->
            roi:    roi
                return [ meta, roi, filtered_callable ]
            no_roi: !roi
                return [ meta, filtered_callable ]
        }

    // Intersect the ROI with the callable regions
    BEDTOOLS_INTERSECT(
        ch_beds_to_intersect.roi,
        ch_fai
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

    def ch_ready_beds = ch_beds_to_intersect.no_roi
        .mix(BEDTOOLS_INTERSECT.out.intersect)

    emit:
    ready_crams = ch_ready_crams    // [ val(meta), path(cram), path(crai) ]
    ready_bams  = ch_ready_bams     // [ val(meta), path(bam), path(bai) ]
    ready_beds  = ch_ready_beds     // [ val(meta), path(bed) ]
    versions    = ch_versions       // [ path(versions) ]
    reports     = ch_reports        // [ path(reports) ]
}
