//
// Concatenate the VCFs back together with bcftools concat
//

include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { TABIX_TABIX     } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_CONCAT_BCFTOOLS {
    take:
        ch_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
        val_tabix // boolean: whether to create a index or not

    main:

    ch_versions = Channel.empty()

    ch_vcfs
        .map { meta, vcf, tbi ->
            def new_meta = meta + [id:meta.sample ?: meta.family]
            [ groupKey(new_meta, meta.split_count), vcf, tbi ]
        }
        .groupTuple()
        .set { ch_concat_input }

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    if(val_tabix) {
        TABIX_TABIX(
            BCFTOOLS_CONCAT.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        BCFTOOLS_CONCAT.out.vcf
            .join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi ->
                // Remove the bed counter from the meta field
                def new_meta = meta - meta.subMap("split_count")
                [ new_meta, vcf, tbi ]
            }
            .set { ch_vcf_tbi }
    } else {
        BCFTOOLS_CONCAT.out.vcf
            .map { meta, vcf ->
                // Remove the bed counter from the meta field
                def new_meta = meta - meta.subMap("split_count")
                [ new_meta, vcf ]
            }
            .set { ch_vcf_tbi }
    }

    emit:
    vcfs = ch_vcf_tbi       // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions  // channel: [ versions.yml ]
}