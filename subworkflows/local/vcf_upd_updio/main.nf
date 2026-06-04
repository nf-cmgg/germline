//
// Run UPDio analysis
//

include { UPDIO             } from '../../../modules/local/updio/main'
include { BCFTOOLS_VIEW     } from '../../../modules/nf-core/bcftools/view'

workflow VCF_UPD_UPDIO {
    take:
        ch_vcfs     // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_peds     // channel: [mandatory] [ val(meta), path(peds) ] => The PED files retrieved from SOMALIER_RELATE
        ch_cnv      // value channel: [optional] [ val(meta), path(cnv) ] => A file with common CNVs to be used by updio
        ch_regions  // value channel: [optional] [ path(bed) ] => A BED file with regions to be used by updio

    main:

    // Filter out all families that are either too small or too big for UPDio analysis
    def ch_trio_vcfs = ch_vcfs
        .filter { meta, _vcf, _tbi ->
            def family_size = meta.family_samples.tokenize(",").size()
            return family_size >= 3 && family_size <= 6
        }

    BCFTOOLS_VIEW(
        ch_trio_vcfs,
        ch_regions,
        [],
        []
    )

    def ch_filter_output = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.index, failOnDuplicate:true, failOnMismatch:true)

    def ch_trio_peds = ch_peds
        .filter { meta, _ped ->
            def family_size = meta.family_samples.tokenize(",").size()
            return family_size >= 3 && family_size <= 6
        }

    def ch_trio_vcfs_family = ch_filter_output
        .join(ch_trio_peds, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi, ped ->
            def meta_list = get_family_data_from_ped(meta, ped)
            [ meta_list, vcf, tbi ]
        }
        .filter { meta, _vcf, _tbi ->
            meta
        }
        .transpose(by:0)

    UPDIO(
        ch_trio_vcfs_family,
        ch_cnv
    )

    emit:
    updio = UPDIO.out.updio    // [ val(meta), path(updio_folder) ]
}

def get_family_data_from_ped(meta, ped) {
    def output = []
    ped.readLines().each { line ->
        if(line.startsWith("#")) { return }
        def split_line = line.split("\t")
        if(split_line[1] != "0" && split_line[2] != "0" && split_line[3] != "0") {
            output.add(meta + [child:split_line[1], father:split_line[2], mother:split_line[3]])
        }
    }
    return output
}
