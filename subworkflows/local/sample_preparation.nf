//
// SAMPLE_PREPARATION
//

include { MERGE_BEDS as MERGE_ROI_PARAMS    } from '../../modules/local/merge_beds'
include { MERGE_BEDS as MERGE_ROI_SAMPLE    } from '../../modules/local/merge_beds'
include { MERGE_BEDS as MERGE_CALLABLE      } from '../../modules/local/merge_beds'
include { SAMTOOLS_MERGE                    } from '../../modules/local/samtools_merge'
include { FILTER_BEDS                       } from '../../modules/local/filter_beds/main'

include { SAMTOOLS_INDEX                    } from '../../modules/nf-core/samtools/index/main'
include { TABIX_TABIX                       } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP as UNZIP_ROI          } from '../../modules/nf-core/tabix/bgzip/main'
include { BEDTOOLS_INTERSECT                } from '../../modules/nf-core/bedtools/intersect/main'
include { MOSDEPTH                          } from '../../modules/nf-core/mosdepth/main'

workflow SAMPLE_PREPARATION {
    take:
        crams             // channel: [mandatory] [ meta, cram, crai ] => sample CRAM files and their optional indices
        roi               // channel: [mandatory] [ meta, roi ] => ROI bed files for WES analysis
        callable          // channel: [mandatory] [ meta, callable ] => callable regions bed files for WES analysis
        fasta             // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai         // channel: [mandatory] [ fasta_fai ] => fasta reference index
        default_roi       // channel: [optional]  [ roi ] => bed containing regions of interest to be used as default

    main:

    ch_versions  = Channel.empty()
    ch_reports   = Channel.empty()

    //
    // Merge the CRAM files if there are multiple per sample
    //

    crams
        .filter { meta, cram, crai ->
            cram != []
        }
        .groupTuple()
        .branch(
            { meta, cram, crai ->
                multiple: cram.size() > 1
                    return [meta, cram]
                single:   cram.size() == 1
                    return [meta, cram, crai]
            }
        )
        .set { cram_branch }

    cram_branch.multiple.dump(tag:'cram_branch_multiple', pretty:true)
    cram_branch.single.dump(tag:'cram_branch_single', pretty:true)

    SAMTOOLS_MERGE(
        cram_branch.multiple,
        fasta,
        fasta_fai
    )

    SAMTOOLS_MERGE.out.cram
        .mix(cram_branch.single
            .map(
                {meta, cram, crai ->
                    [ meta, cram[0], crai[0]]
                }
            )
        )
        .branch(
            { meta, cram, crai=[] ->
                not_indexed: crai == []
                    return [ meta, cram ]
                indexed: crai != []
                    return [ meta, cram, crai ]
            }
        )
        .set { merged_crams }

    merged_crams.not_indexed
        .tap { crams_to_index }
        .dump(tag:'merged_crams_not_indexed', pretty:true)
        .set { crams_without_index }
    merged_crams.indexed.dump(tag:'merged_crams_indexed', pretty:true)

    SAMTOOLS_INDEX(
        crams_to_index
    )

    crams_without_index
        .join(SAMTOOLS_INDEX.out.crai, failOnDuplicate: true, failOnMismatch: true)
        .mix(merged_crams.indexed)
        .tap { mosdepth_crams }
        .dump(tag:'ready_crams', pretty:true)
        .set { ready_crams }

    //
    // Preprocess the ROI BED files => merge overlapping 
    //

    roi
        .groupTuple() // A specified size isn't needed here since this runs before any process using ROI files is executed
        .branch { meta, roi ->
            // Determine whether or not there is an ROI file given to the current sample
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
        .set { roi_branch }

    // Merge the ROI BED files if multiple samples are given, also merges overlapping regions in the BED files
    MERGE_ROI_SAMPLE(
        roi_branch.found
    )
    ch_versions = ch_versions.mix(MERGE_ROI_SAMPLE.out.versions.first())

    // Add the default ROI file to all samples without an ROI file 
    // if an ROI BED file has been given through the --roi parameter
    if (default_roi) {
        MERGE_ROI_PARAMS(
            default_roi.map { [[id:"default_roi"], it]}
        )
        ch_versions = ch_versions.mix(MERGE_ROI_PARAMS.out.versions)

        roi_branch.missing
            .groupTuple() // A specified size isn't needed here since this runs before any process using the default ROI file is executed
            .combine(MERGE_ROI_PARAMS.out.bed.map { it[1] })
            .map { meta, missing, default_roi ->
                [ meta, default_roi ]
            }
            .set { missing_rois }
    } else {
        roi_branch.missing.set { missing_rois }
    }

    missing_rois
        .mix(MERGE_ROI_SAMPLE.out.bed)
        .set { ready_rois }

    // Merge callable BED files
    callable
        .groupTuple() // A specified size isn't needed here since this runs before any process using the callable BED files is executed
        .branch { meta, callable ->
            // Determine whether or not there is a callable file given to the current sample
            def is_present = false
            def output_callable = []
            for( entry : callable ) {
                if(entry != []){
                    output_callable.add(entry)
                    is_present = true
                }
            }
            found:      output_callable.size() >= 2
                return [ meta, output_callable ]
            missing:    output_callable.size() < 2
                return [ meta, [] ]
        }
        .set { callable_branch }

    MERGE_CALLABLE(
        callable_branch.found
    )
    ch_versions = ch_versions.mix(MERGE_CALLABLE.out.versions)

    callable_branch.missing
        .mix(MERGE_CALLABLE.out.bed)
        .set { ready_callable }

    // Mix the default ROI with the sample-specific ROI files
    // and determine whether the sample is WES or WGS (aka whether they 
    // have an ROI file or not)
    missing_rois
        .mix(MERGE_ROI_SAMPLE.out.bed)
        .join(ready_callable, failOnMismatch:true, failOnDuplicate:true)
        .set { ready_regions }

    //
    // Create callable regions if they don't exist
    //

    // Check if a BED with callable regions has been given for each sample
    ready_regions
        .branch { meta, roi, callable ->
            present:        callable
            not_present:    !callable
        }
        .set { callable_branch }

    callable_branch.present.view()

    // Create BEDs with callable regions using Mosdepth
    callable_branch.not_present
        .join(ready_crams, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, roi, callable, cram, crai ->
            [ meta, cram, crai, roi ]
        }
        .set { mosdepth_input }

    MOSDEPTH(
        mosdepth_input,
        fasta.map { [[], it] }
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    // Join all channels back together
    MOSDEPTH.out.quantized_bed
        .join(callable_branch.not_present, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, new_callable, roi, callable ->
            [ meta, roi, new_callable ]
        }
        .mix(callable_branch.present)
        .set { beds_to_filter }

    // Filter out the regions with no coverage
    FILTER_BEDS(
        beds_to_filter.map { meta, roi, callable -> [ meta, callable ]}
    )
    ch_versions = ch_versions.mix(FILTER_BEDS.out.versions)

    FILTER_BEDS.out.bed
        .join(beds_to_filter, failOnDuplicate:true, failOnMismatch:true)
        .branch { meta, filtered_callable, roi, callable ->
            roi:    roi
                return [ meta, roi, filtered_callable]
            no_roi: !roi
                return [ meta, filtered_callable ]
        }
        .set { beds_to_intersect }

    // Intersect the ROI with the callable regions
    BEDTOOLS_INTERSECT(
        beds_to_intersect.roi,
        fasta_fai.map { [[], it] }
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

    beds_to_intersect.no_roi
        .mix(BEDTOOLS_INTERSECT.out.intersect)
        .set { ready_beds }

    emit:
    ready_crams
    ready_beds
    versions = ch_versions
    reports  = ch_reports
}
