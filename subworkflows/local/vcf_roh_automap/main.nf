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

    def hg_genome = val_genome == "GRCh38" ? "hg38" : val_genome == "GRCh37" ? "hg19" : val_genome

    if(!["hg38", "hg19"].contains(hg_genome)) {
        error("Genome '${hg_genome}' is not supported by automap. Available options are hg38/GRCh38 or hg19/GRCh37.")
    }

    // Merge the repeat BED files from the container if no container has been given
    def ch_valid_repeats = channel.empty()
    if (!ch_repeats) {
        AUTOMAP_REPEATS(
            channel.value([[id:"${val_genome}_repeats"], val_genome])
        )
        ch_valid_repeats = AUTOMAP_REPEATS.out.repeats.collect()
    } else {
        ch_valid_repeats = ch_repeats
    }

    def ch_automap_input = ch_vcfs.filter { meta, vcf, tbi ->
        // Check the amount of variants for all VCFs that are under 1 MB (automap requires at least 10000 variants to run)
        if(vcf?.size() < 1000000) {
            def var_count = 0
            vcf.withInputStream { stream ->
                new java.util.zip.GZIPInputStream(stream).withReader('UTF-8') { reader ->
                    reader.eachLine { line ->
                        if(!line.startsWith("#")) {
                            var_count += 1
                        }
                    }
                }
            }
            if(var_count < 10000) {
                log.warn("VCF file '${vcf}' for sample '${meta.id}' only has ${var_count} variants. Automap requires at least 10000 variants to run. Skipping this VCF.")
                return false
            }
        }
        return true
    }

    AUTOMAP_AUTOMAP(
        ch_automap_input,
        ch_valid_repeats,
        ch_panel,
        hg_genome
    )

    emit:
    automap = AUTOMAP_AUTOMAP.out.automap   // [ val(meta), path(automap_folder) ]
}
