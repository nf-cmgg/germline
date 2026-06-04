//
// Filter the VCFs
//

include { BCFTOOLS_FILTER } from '../../../modules/nf-core/bcftools/filter/main'

workflow VCF_FILTER_BCFTOOLS {
    take:
        ch_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    main:

    BCFTOOLS_FILTER(
        ch_vcfs
    )

    def ch_filter_vcfs = BCFTOOLS_FILTER.out.vcf
        .join(BCFTOOLS_FILTER.out.index, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vcfs = ch_filter_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]
}
