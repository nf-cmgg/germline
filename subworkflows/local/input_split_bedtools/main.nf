//
// Split the BED files and return a channel ready to use for scatter/gather flows
//

include { BEDTOOLS_SPLIT } from '../../../modules/nf-core/bedtools/split/main'

workflow INPUT_SPLIT_BEDTOOLS {
    take:
        ch_beds     // channel: [mandatory] [ val(meta), path(bed) ] => bed files
        ch_inputs   // channel: [mandatory] [ val(meta), path(input), path(input_index), path(extra_file) ] => input files

    main:

    ch_versions = Channel.empty()

    BEDTOOLS_SPLIT(
        ch_beds.map { it + [params.scatter_count]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    ch_inputs
        .join(BEDTOOLS_SPLIT.out.beds, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, input, input_index, extra_file, beds ->
            // Determine the amount of BED files per sample
            bed_is_list = beds instanceof ArrayList
            new_meta = meta + [split_count: bed_is_list ? beds.size() : 1]
            [ new_meta, input, input_index, extra_file, bed_is_list ? beds : [beds] ]
        }
        .transpose(by:3) // Create one channel entry for each BED file per sample
        .map { meta, input, input_index, extra_file, bed ->
            // Set the base name of the BED file as the ID (this will look like sample_id.xxxx, where xxxx are numbers)
            new_meta = meta + [id:bed.baseName]
            [ new_meta, input, input_index, extra_file, bed ]
        }
        .set { ch_split_output }

    emit:
        split = ch_split_output // channel: [ val(meta), path(input), path(input_index), path(extra_file), path(bed) ]
        versions = ch_versions  // channel: [ versions.yml ]
}