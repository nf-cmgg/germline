include { SCATTER_BEDS  } from '../../modules/local/scatter_beds'

workflow BED_SCATTER_GROOVY {

    take:
    ch_regions      // channel: [ meta, regions_file ]
    minimum_size    // integer: The minimum size each BED file should contain

    main:

    ch_versions = Channel.empty()

    SCATTER_BEDS(
        ch_regions,
        minimum_size
    )

    SCATTER_BEDS.out.scatter
        .map { meta, beds ->
            [ meta, beds, beds instanceof ArrayList ? beds.size() : 1 ]
        }
        .transpose()
        .dump(tag:'scatter_beds', pretty:true)
        .set { ch_regions_files }

    emit:
    scattered = ch_regions_files    // channel: [ meta, scattered_region_file, scatter_count]

    versions = ch_versions          // channel: [ versions.yml ]
}

