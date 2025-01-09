//
// Run UPDio analysis
//

include { AUTOMAP_REPEATS    } from '../../../modules/local/automap/repeats/main'
include { AUTOMAP_AUTOMAP    } from '../../../modules/local/automap/automap/main'

workflow VCF_ROH_AUTOMAP {
    take:
        ch_vcfs     // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_repeats  // value channel: [optional] [ val(meta), path(repeats) ] => A BED file with repeat regions in the genome
        ch_panel    // value channel: [optional]  [ val(meta), path(panel) ] => A TXT file containing the locations of the genes to be analyzed
        val_genome  // value: [mandatory] => The genome to be used by automap
    main:

    def ch_versions = Channel.empty()

    def hg_genome = val_genome == "GRCh38" ? "hg38" : val_genome == "GRCh37" ? "hg19" : val_genome

    if(!["hg38", "hg19"].contains(hg_genome)) {
        error("Genome '${hg_genome}' is not supported by automap. Available options are hg38/GRCh38 or hg19/GRCh37.")
    }

    // Merge the repeat BED files from the container if no container has been given
    def ch_valid_repeats = Channel.empty()
    if (!ch_repeats) {
        AUTOMAP_REPEATS(
            Channel.value([[id:"${val_genome}_repeats"], val_genome])
        )
        ch_versions = ch_versions.mix(AUTOMAP_REPEATS.out.versions)
        ch_valid_repeats = AUTOMAP_REPEATS.out.repeats.collect()
    } else {
        ch_valid_repeats = ch_repeats
    }

    AUTOMAP_AUTOMAP(
        ch_vcfs,
        ch_valid_repeats,
        ch_panel,
        hg_genome
    )
    ch_versions = ch_versions.mix(AUTOMAP_AUTOMAP.out.versions.first())

    emit:
    automap = AUTOMAP_AUTOMAP.out.automap   // [ val(meta), path(automap_folder) ]
    versions = ch_versions                  // [ path(versions) ]
}
