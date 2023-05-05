//
// ADD PED HEADER
//

include { RTGTOOLS_PEDFILTER    } from '../../modules/nf-core/rtgtools/pedfilter/main'
include { BCFTOOLS_ANNOTATE     } from '../../modules/nf-core/bcftools/annotate/main'

workflow ADD_PED_HEADER {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf) ] => The post-processed VCFs
        ch_somalier_samples_tsv // channel: [mandatory] [ val(meta), path(samples_tsv) ] => The samples TSV retrieved from SOMALIER_RELATE

    main:

    ch_versions = Channel.empty()

    //
    // Remove extra columns from the samples TSV and convert to a VCF header
    //

    RTGTOOLS_PEDFILTER(
        ch_somalier_samples_tsv
    )
    ch_versions = ch_versions.mix(RTGTOOLS_PEDFILTER.out.versions.first())

    //
    // Add the PED headers to the VCF using bcftools annotate --header-lines
    //

    ch_vcfs
        .join(RTGTOOLS_PEDFILTER.out.output, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, ped_vcf ->
            [ meta, vcf, [], [], [], ped_vcf ]
        }
        .set { ch_annotate_input }
    

    BCFTOOLS_ANNOTATE(
        ch_annotate_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    emit:
    ped_vcfs = BCFTOOLS_ANNOTATE.out.vcf // [ val(meta), path(vcf) ]
    versions = ch_versions               // [ path(versions) ]
}
