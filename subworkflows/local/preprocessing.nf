//
// PREPROCESSING
//

include { MERGE_BEDS                } from '../../modules/local/merge_beds'
include { SAMTOOLS_MERGE            } from '../../modules/local/samtools_merge'

include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'
include { MOSDEPTH                  } from '../../modules/nf-core/mosdepth/main'
include { TABIX_BGZIP as UNZIP_BEDS } from '../../modules/nf-core/tabix/bgzip/main'

workflow PREPROCESSING {
    take:
        crams             // channel: [mandatory] [ meta, cram, crai ] => sample CRAM files and their optional indices
        beds              // channel: [mandatory] [ meta, bed ] => bed files
        fasta             // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai         // channel: [mandatory] [ fasta_fai ] => fasta reference index

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
        .tap { mosdept_crams }
        .dump(tag:'ready_crams', pretty:true)
        .set { ready_crams }

    //
    // Merge the BED files if there are multiple per sample
    //

    beds
        .groupTuple()
        .branch(
            { meta, bed ->
                multiple: bed.size() > 1
                    return [meta, bed]
                single:   bed.size() == 1 && bed != [[]]
                    return [meta, bed[0]]
                no_bed: bed == [[]]
                    return [meta, []]
            }
        )
        .set { bed_branch }

    MERGE_BEDS(
        bed_branch.multiple
    )

    bed_branch.no_bed
        .join(mosdept_crams)
        .map { meta, bed, cram, crai ->
            [ meta, cram, crai ]
        }
        .set { mosdepth_input }

    //
    // Create BED files if none are supplied
    //

    MOSDEPTH(
        mosdepth_input,
        [[],[]],
        fasta.map { [[],it]}
    )

    UNZIP_BEDS(
        MOSDEPTH.out.quantized_bed
    )

    MERGE_BEDS.out.bed
        .mix(UNZIP_BEDS.out.output)
        .mix(bed_branch.single)
        .dump(tag:'merged_beds_preprocessing', pretty:true)
        .set { ready_beds }


    emit:
    ready_crams
    ready_beds
    versions = ch_versions
    reports  = ch_reports
}
