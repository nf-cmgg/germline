//
// ADD PED HEADER
//

include { RTGTOOLS_PEDFILTER    } from '../../modules/nf-core/rtgtools/pedfilter/main'
include { BCFTOOLS_REHEADER     } from '../../modules/nf-core/bcftools/reheader/main'

include { MERGE_HEADERS         } from '../../modules/local/merge_headers'

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
    // Create the new headers and replace the VCF headers with the new headers
    //

    // TODO: rewrite this module in a better language than groovy
    MERGE_HEADERS(
        ch_vcfs.join(RTGTOOLS_PEDFILTER.out.output, failOnDuplicate: true, failOnMismatch: true)
    )
    ch_versions = ch_versions.mix(MERGE_HEADERS.out.versions.first())

    ch_vcfs
        .join(MERGE_HEADERS.out.header, failOnDuplicate: true, failOnMismatch: true)
        .set { ch_bcftools_reheader_input }

    BCFTOOLS_REHEADER(
        ch_bcftools_reheader_input,
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions.first())

    emit:
    ped_vcfs = BCFTOOLS_REHEADER.out.vcf // [ val(meta), path(vcf) ]
    versions = ch_versions               // [ path(versions) ]
}
