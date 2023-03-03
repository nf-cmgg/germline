//
// ADD PED HEADER
//

include { RTGTOOLS_PEDFILTER    } from '../../modules/nf-core/rtgtools/pedfilter/main'
include { BCFTOOLS_REHEADER     } from '../../modules/nf-core/bcftools/reheader/main'

include { MERGE_HEADERS         } from '../../modules/local/merge_headers'

workflow ADD_PED_HEADER {
    take:
        vcfs                 // channel: [mandatory] [ meta, vcfs ] => The post-processed VCFs
        somalier_samples_tsv // channel: [mandatory] [ meta, samples_tsv ] => The samples TSV retrieved from SOMALIER_RELATE

    main:

    ch_versions         = Channel.empty()

    //
    // Remove extra columns from the samples TSV and convert to a VCF header
    //

    RTGTOOLS_PEDFILTER(
        somalier_samples_tsv
    )

    ch_versions = ch_versions.mix(RTGTOOLS_PEDFILTER.out.versions)

    //
    // Create the new headers
    //

    MERGE_HEADERS(
        vcfs.join(RTGTOOLS_PEDFILTER.out.output, failOnDuplicate: true, failOnMismatch: true)
    )

    ch_versions = ch_versions.mix(MERGE_HEADERS.out.versions)

    vcfs
        .join(MERGE_HEADERS.out.header, failOnDuplicate: true, failOnMismatch: true)
        .set { bcftools_reheader_input }

    BCFTOOLS_REHEADER(
        bcftools_reheader_input,
        []
    )

    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    emit:
    ped_vcfs = BCFTOOLS_REHEADER.out.vcf
    versions = ch_versions
}
