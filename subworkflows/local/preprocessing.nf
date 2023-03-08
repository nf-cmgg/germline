//
// PREPROCESSING
//

include { MERGE_BEDS as MERGE_ROI       } from '../../modules/local/merge_beds'
include { SAMTOOLS_MERGE                } from '../../modules/local/samtools_merge'

include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { GOLEFT_INDEXSPLIT             } from '../../modules/nf-core/goleft/indexsplit/main'
include { TABIX_BGZIP as UNZIP_BEDS     } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as UNZIP_ROI      } from '../../modules/nf-core/tabix/bgzip/main'
include { BEDTOOLS_INTERSECT            } from '../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_MERGE                } from '../../modules/nf-core/bedtools/merge/main'

workflow PREPROCESSING {
    take:
        crams             // channel: [mandatory] [ meta, cram, crai ] => sample CRAM files and their optional indices
        roi               // channel: [mandatory] [ meta, bed ] => bed files containing regions of interest
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
        .join(SAMTOOLS_INDEX.out.crai)
        .mix(merged_crams.indexed)
        .tap { mosdepth_crams }
        .dump(tag:'ready_crams', pretty:true)
        .set { ready_crams }

    //
    // Unzip the Gzipped beds
    //

    roi
        .branch { meta, bed ->
            gunzipped: bed == [] || bed.getExtension() == "bed"
            gzipped: bed.getExtension() == "gz"
        }
        .set { roi_branch }

    UNZIP_ROI(
        roi_branch.gzipped
    )

    ch_versions = ch_versions.mix(UNZIP_ROI.out.versions)

    //
    // Merge the BED files if there are multiple per sample
    //

    roi_branch.gunzipped
        .mix(UNZIP_ROI.out.output)
        .groupTuple() //TODO add a size here
        .branch(
            { meta, bed ->
                multiple: bed.size() > 1
                    return [meta, bed]
                single:   bed.size() == 1
                    return [meta, bed[0]]
            }
        )
        .set { merge_roi_branch }

    MERGE_ROI(
        merge_roi_branch.multiple
    )

    ch_versions = ch_versions.mix(MERGE_ROI.out.versions)

    merge_roi_branch.single
        .mix(MERGE_ROI.out.bed)
        .set { ready_roi }

    //
    // Create BED files with bins containing equal amounts of reads
    //

    ready_crams
        .map { meta, cram, crai ->
            [ meta, crai ]
        }
        .set { indexsplit_input }

    GOLEFT_INDEXSPLIT(
        indexsplit_input,
        fasta_fai.map { [[], it] },
        params.scatter_count
    )

    ch_versions = ch_versions.mix(GOLEFT_INDEXSPLIT.out.versions)

    //
    // Intersect the ROI BEDs and binned region BEDs
    //

    ready_roi
        .join(GOLEFT_INDEXSPLIT.out.bed)
        .combine(default_roi)
        .branch { meta, roi, bed, default_roi ->
            out_roi = roi ?: default_roi ?: []
            intersect: out_roi != []
                return [ meta, out_roi, bed ]
            no_intersect: out_roi == []
                return [ meta, bed ]
        }
        .set { intersect_branch }

    BEDTOOLS_INTERSECT(
        intersect_branch.intersect,
        "bed"
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

    BEDTOOLS_MERGE(
        BEDTOOLS_INTERSECT.out.intersect
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

    intersect_branch.no_intersect
        .mix(BEDTOOLS_MERGE.out.bed)
        .set { ready_beds }

    emit:
    ready_crams
    ready_beds
    versions = ch_versions
    reports  = ch_reports
}
