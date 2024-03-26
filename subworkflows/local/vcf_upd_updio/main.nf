//
// Run UPDio analysis
//

include { UPDIO    } from '../../../modules/local/updio/main'

workflow VCF_UPD_UPDIO {
    take:
        ch_vcfs // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_peds // channel: [mandatory] [ val(meta), path(peds) ] => The PED files retrieved from SOMALIER_RELATE
        ch_cnv  // value channel: [optional] [ val(meta), path(cnv) ] => A file with common CNVs to be used by updio

    main:

    ch_versions = Channel.empty()

    // Filter out all families that don't have exactly 3 samples
    ch_vcfs
        .filter { it[0].family_count == 3 }
        .set { ch_trio_vcfs }

    ch_peds
        .filter { it[0].family_count == 3 }
        .set { ch_trio_peds }

    ch_trio_vcfs
        .join(ch_trio_peds)
        .map { meta, vcf, tbi, ped ->
            def new_meta = get_family_data_from_ped(meta, ped)
            [ new_meta, vcf, tbi ]
        }
        .filter { it[0].child != null && it[0].mother != null && it[0].father != null }
        .set { ch_trio_vcfs_family }

    UPDIO(
        ch_trio_vcfs_family,
        ch_cnv
    )

    emit:
    // ped_vcfs = BCFTOOLS_ANNOTATE.out.vcf    // [ val(meta), path(vcf) ]
    versions = ch_versions                  // [ path(versions) ]
}

def get_family_data_from_ped(meta, ped) {
    def child = null
    def mother = null
    def father = null
    ped.readLines().each { line ->
        if(line.startsWith("#")) { return }
        def split_line = line.split("\t")
        if(split_line[1] != "0" && split_line[2] != "0" && split_line[3] != "0") {
            child = split_line[1]
            father = split_line[2]
            mother = split_line[3]
        }
    }
    return meta + [child:child, mother:mother, father:father]
}
