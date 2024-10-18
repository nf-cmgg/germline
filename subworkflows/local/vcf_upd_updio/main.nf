//
// Run UPDio analysis
//

include { UPDIO             } from '../../../modules/local/updio/main'
include { BCFTOOLS_FILTER   } from '../../../modules/nf-core/bcftools/filter'

workflow VCF_UPD_UPDIO {
    take:
        ch_vcfs // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_peds // channel: [mandatory] [ val(meta), path(peds) ] => The PED files retrieved from SOMALIER_RELATE
        ch_cnv  // value channel: [optional] [ val(meta), path(cnv) ] => A file with common CNVs to be used by updio

    main:

    def ch_versions = Channel.empty()

    // Filter out all families that have less than 3 samples
    def ch_trio_vcfs = ch_vcfs
        .filter { meta, _vcf, _tbi ->
            meta.family_samples.tokenize(",").size() >= 3
        }

    BCFTOOLS_FILTER(
        ch_trio_vcfs
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    def ch_filter_output = BCFTOOLS_FILTER.out.vcf
        .join(BCFTOOLS_FILTER.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    def ch_trio_peds = ch_peds
        .filter { meta, _ped ->
            meta.family_samples.tokenize(",").size() >= 3
        }

    def ch_trio_vcfs_family = CustomChannelOperators.joinOnKeys(
            [failOnDuplicate:true, failOnMismatch:true],
            ch_filter_output,
            ch_trio_peds,
            ["id", "family", "family_samples", "caller"]
        )
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
    ch_versions = ch_versions.mix(UPDIO.out.versions.first())

    emit:
    updio = UPDIO.out.updio    // [ val(meta), path(updio_folder) ]
    versions = ch_versions     // [ path(versions) ]
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
