//
// ADD PED HEADER
//

include { RTGTOOLS_PEDFILTER    } from '../../../modules/nf-core/rtgtools/pedfilter/main'
include { BCFTOOLS_ANNOTATE     } from '../../../modules/nf-core/bcftools/annotate/main'

workflow VCF_PED_RTGTOOLS {
    take:
        ch_vcfs // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_peds // channel: [mandatory] [ val(meta), path(peds) ] => The PED files retrieved from SOMALIER_RELATE

    main:

    //
    // Remove extra columns from the samples TSV and convert to a VCF header
    //

    RTGTOOLS_PEDFILTER(
        ch_peds
    )

    //
    // Add the PED headers to the VCF using bcftools annotate --header-lines
    //

    def ch_annotate_input = ch_vcfs
        .join(RTGTOOLS_PEDFILTER.out.output, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi, ped_vcf ->
            [ meta, vcf, tbi, [], [], [], ped_vcf, [] ]
        }

    BCFTOOLS_ANNOTATE(
        ch_annotate_input
    )

    def ch_ped_vcfs = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.index, failOnDuplicate:true, failOnMismatch:true)

    emit:
    ped_vcfs = ch_ped_vcfs  // [ val(meta), path(vcf), path(tbi) ]
}
