//
// SAMPLE_PREPARATION
//

include { MERGE_BEDS as MERGE_ROI_PARAMS    } from '../../modules/local/merge_beds'
include { MERGE_BEDS as MERGE_ROI_SAMPLE    } from '../../modules/local/merge_beds'
include { FILTER_BEDS                       } from '../../modules/local/filter_beds/main'

include { SAMTOOLS_MERGE                    } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                    } from '../../modules/nf-core/samtools/index/main'
include { TABIX_TABIX                       } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP as UNZIP_ROI          } from '../../modules/nf-core/tabix/bgzip/main'
include { BEDTOOLS_INTERSECT                } from '../../modules/nf-core/bedtools/intersect/main'
include { MOSDEPTH                          } from '../../modules/nf-core/mosdepth/main'

workflow SAMPLE_PREPARATION {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their optional indices
        ch_roi               // channel: [optional]  [ val(meta), path(roi) ] => ROI bed files for WES analysis
        ch_fasta             // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_default_roi       // channel: [optional]  [ path(roi) ] => bed containing regions of interest to be used as default

    main:

    ch_versions  = Channel.empty()
    ch_reports   = Channel.empty()

    //
    // Merge the CRAM files if there are multiple per sample
    //

    ch_crams
        .groupTuple() // No size needed here because this runs before any process
        .branch(
            { meta, cram, crai ->
                multiple: cram.size() > 1
                    return [meta, cram]
                single:   cram.size() == 1
                    return [meta, cram[0], crai[0]]
            }
        )
        .set { ch_cram_branch }

    ch_cram_branch.multiple.dump(tag:'cram_branch_multiple', pretty:true)
    ch_cram_branch.single.dump(tag:'cram_branch_single', pretty:true)

    SAMTOOLS_MERGE(
        ch_cram_branch.multiple,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    //
    // Index the CRAM files which have no index
    //

    SAMTOOLS_MERGE.out.cram
        .mix(ch_cram_branch.single)
        .branch { meta, cram, crai=[] ->
            not_indexed: crai == []
                return [ meta, cram ]
            indexed: crai != []
                return [ meta, cram, crai ]
        }
        .set { ch_merged_crams }

    ch_merged_crams.not_indexed.dump(tag:'merged_crams_not_indexed', pretty:true)
    ch_merged_crams.indexed.dump(tag:'merged_crams_indexed', pretty:true)

    SAMTOOLS_INDEX(
        ch_merged_crams.not_indexed
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_merged_crams.not_indexed
        .join(SAMTOOLS_INDEX.out.crai, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_merged_crams.indexed)
        .dump(tag:'ready_crams', pretty:true)
        .set { ch_ready_crams }

    //
    // Preprocess the ROI BED files => sort and merge overlapping regions
    //

    ch_roi
        .groupTuple() // A specified size isn't needed here since this runs before any process using ROI files is executed
        .branch { meta, roi ->
            // Determine whether there is an ROI file given to the current sample
            // It's possible that a sample is given multiple times in the samplesheet, in which
            // case they have been merged earlier. This code checks if at least one entry of the same
            // sample contains an ROI file
            def is_present = false
            def output_roi = []
            for( entry : roi) {
                if(entry != []){
                    output_roi.add(entry)
                    is_present = true
                }
            }
            found:      is_present
                return [ meta, output_roi ]
            missing:    !is_present
                return [ meta, [] ]
        }
        .set { ch_roi_branch }

    MERGE_ROI_SAMPLE(
        ch_roi_branch.found,
        ch_fai
    )
    ch_versions = ch_versions.mix(MERGE_ROI_SAMPLE.out.versions.first())

    // Add the default ROI file to all samples without an ROI file 
    // if an ROI BED file has been given through the --roi parameter
    if (ch_default_roi) {
        MERGE_ROI_PARAMS(
            ch_default_roi.map { [[id:"default_roi"], it]},
            ch_fai
        )
        ch_versions = ch_versions.mix(MERGE_ROI_PARAMS.out.versions)

        ch_roi_branch.missing
            .groupTuple() // A specified size isn't needed here since this runs before any process using the default ROI file is executed
            .combine(MERGE_ROI_PARAMS.out.bed.map { it[1] })
            .map { meta, missing, default_roi ->
                [ meta, default_roi ]
            }
            .set { ch_missing_rois }
    } else {
        ch_roi_branch.missing.set { ch_missing_rois }
    }

    ch_missing_rois
        .mix(MERGE_ROI_SAMPLE.out.bed)
        .set { ch_ready_rois }

    //
    // Create callable regions
    //

    ch_ready_crams
        .join(ch_ready_rois, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_mosdepth_input }

    MOSDEPTH(
        ch_mosdepth_input,
        ch_fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    ch_ready_rois
        .join(MOSDEPTH.out.quantized_bed, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_beds_to_filter }

    // Filter out the regions with no coverage
    FILTER_BEDS(
        ch_beds_to_filter.map { meta, roi, callable -> [ meta, callable ]}
    )
    ch_versions = ch_versions.mix(FILTER_BEDS.out.versions)

    FILTER_BEDS.out.bed
        .join(ch_beds_to_filter, failOnDuplicate:true, failOnMismatch:true)
        .branch { meta, filtered_callable, roi, callable ->
            roi:    roi
                return [ meta, roi, filtered_callable ]
            no_roi: !roi
                return [ meta, filtered_callable ]
        }
        .set { ch_beds_to_intersect }

    // Intersect the ROI with the callable regions
    BEDTOOLS_INTERSECT(
        ch_beds_to_intersect.roi,
        ch_fai
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

    ch_beds_to_intersect.no_roi
        .mix(BEDTOOLS_INTERSECT.out.intersect)
        .set { ch_ready_beds }

    emit:
    ready_crams = ch_ready_crams    // [ val(meta), path(cram), path(crai) ]
    ready_beds  = ch_ready_beds     // [ val(meta), path(bed) ]
    versions    = ch_versions       // [ path(versions) ]
    reports     = ch_reports        // [ path(reports) ]
}
