//
// Filter the VCFs
//

include { BCFTOOLS_FILTER as FILTER_1 } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_2 } from '../../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX                 } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_FILTER_BCFTOOLS {
    take:
        ch_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
        val_tabix // boolean: whether to create a index or not

    main:

    def ch_versions = Channel.empty()

    FILTER_1(
        ch_vcfs.map { meta, vcf, tbi=[] -> [ meta, vcf, tbi ]}
    )
    ch_versions = ch_versions.mix(FILTER_1.out.versions.first())

    FILTER_2(
        FILTER_1.out.vcf.map { meta, vcf -> [ meta, vcf, [] ]}
    )
    ch_versions = ch_versions.mix(FILTER_2.out.versions.first())

    def ch_filter_vcfs = Channel.empty()
    if(val_tabix) {
        TABIX_TABIX(
            FILTER_2.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        ch_filter_vcfs = FILTER_2.out.vcf
            .join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
    } else {
        ch_filter_vcfs = FILTER_2.out.vcf
    }


    emit:
    vcfs = ch_filter_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions // channel: [ versions.yml ]
}
