//
// VCF_QC
//

include { BCFTOOLS_STATS                  } from '../../modules/nf-core/modules/bcftools/stats/main'
include { VCFTOOLS as VCFTOOLS_SUMMARY    } from '../../modules/nf-core/modules/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_COUNT } from '../../modules/nf-core/modules/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_QUAL  } from '../../modules/nf-core/modules/vcftools/main'

workflow VCF_QC {
    take:
        vcfs   // channel: [mandatory] vcfs

    main:

    ch_versions         = Channel.empty()

    //
    // Perform all quality control steps
    //

    BCFTOOLS_STATS(
        vcfs
    )
    VCFTOOLS_TSTV_COUNT(
        vcfs, 
        [], 
        []
    )
    VCFTOOLS_TSTV_QUAL(
        vcfs, 
        [], 
        []
    )
    VCFTOOLS_SUMMARY(
        vcfs, 
        [], 
        []
    )

    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_TSTV_COUNT.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_TSTV_QUAL.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_SUMMARY.out.versions)

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_tstv_count     = VCFTOOLS_TSTV_COUNT.out.tstv_count
    vcftools_tstv_qual      = VCFTOOLS_TSTV_QUAL.out.tstv_qual
    vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary
    versions                = ch_versions
}