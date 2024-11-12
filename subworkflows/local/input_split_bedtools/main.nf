//
// Split the BED files and return a channel ready to use for scatter/gather flows
//

include { BEDTOOLS_SPLIT } from '../../../modules/nf-core/bedtools/split/main'

workflow INPUT_SPLIT_BEDTOOLS {
    take:
        ch_beds     // channel: [mandatory] [ val(meta), path(bed), val(scatter_count) ] => bed files
        ch_inputs   // channel: [mandatory] [ val(meta), path(input), path(input_index) ] => input files

    main:

    def ch_versions = Channel.empty()

    BEDTOOLS_SPLIT(
        ch_beds
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    def ch_split_output = ch_inputs
        .join(BEDTOOLS_SPLIT.out.beds, failOnDuplicate: true, failOnMismatch: true)
        .map { row ->
            def meta = row[0]
            def beds = row[-1]
            // Determine the amount of BED files per sample
            def bed_is_list = beds instanceof ArrayList
            def new_meta = meta + [split_count: bed_is_list ? beds.size() : 1]
            def bed_output = bed_is_list ? [beds] : [[beds]]
            return [new_meta] + bed_output + row[1..-2]
        }
        .transpose(by:1) // Create one channel entry for each BED file per sample
        .map { row ->
            // Set the base name of the BED file as the ID (this will look like sample_id.xxxx, where xxxx are numbers)
            def new_meta = row[0] + [id:row[1].baseName]
            return [ new_meta ] + row[2..-1] + [ row[1] ]
        }

    emit:
    split = ch_split_output // channel: [ val(meta), path(input), path(input_index), path(bed) ]
    versions = ch_versions  // channel: [ versions.yml ]
}
