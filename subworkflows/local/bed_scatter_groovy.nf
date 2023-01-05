workflow BED_SCATTER_GROOVY {

    take:
    ch_regions      // channel: [ meta, regions_file ]
    minimum_size    // integer: The minimum size each BED file should contain

    main:

    ch_versions = Channel.empty()

    ch_regions
        .multiMap { meta, regions_file ->
            file_type = regions_file.getExtension()
            scattered_regions = scatter_regions(meta, regions_file, minimum_size)
            total = scattered_regions.size()
            range = 1..total

            regions: [ meta, scattered_regions, range, file_type ]
            meta:    [ meta.id, meta, total ]
        }
    .set { ch_scattered }

    ch_scattered.regions
        .dump(tag:'scatter_regions', pretty:true)
        .transpose()
        .collectFile()
        { meta, regions, number, file_type ->
            [ "${meta.id}_${number}.${file_type}", regions.join("\n") ]
        }
        .map { file ->
            id = file.baseName.tokenize("_")[0..-2].join("_")
            [ id, file ]
        }
        .dump(tag:'scatter_beds', pretty:true)
        .combine(
            ch_scattered.meta.dump(tag:'scatter_meta', pretty:true),
            by:0
        )
        .map { id, file, meta, total ->
            [ meta, file, total ]
        }
        .dump(tag:'scatter_output', pretty:true)
        .set { ch_regions_files }

    emit:
    scattered = ch_regions_files    // channel: [ meta, scattered_region_file, scatter_count]

    versions = ch_versions          // channel: [ versions.yml ]
}

def scatter_regions(meta, file, minimum_size) {
    stringified_file = file.toString().toLowerCase()
    output_regions = []
    if (stringified_file.endsWith("bed")) {
        file.withReader {
            String line
            String status = "Added"
            Integer file_count = 0
            Integer current_size = 0
            current_regions = []
            while( line = it.readLine() ) {
                if( line ==~ /^#.*$/) {continue} // skip comments
                status = "Not in there"
                entry = line.tokenize("\t")
                current_regions.add(line)
                current_size = current_size + (entry[2].toInteger() - entry[1].toInteger())

                if( current_size >= minimum_size ){
                    output_regions.add(current_regions)

                    status = "Added"
                    current_size = 0
                    current_regions = []
                }
            }
            if (status == "Not in there") {
                output_regions.add(current_regions)
            }
        }
    }

    // Check for stubs!
    if(!output_regions){
        println("WARNING: No BED contents detected when scattering, created some artificial regions for ${meta.id} - ignore this when using stub runs")
        output_regions = [["chr20\t100\t200"],["chr20\t300\t400"]]
    }

    return output_regions
}
