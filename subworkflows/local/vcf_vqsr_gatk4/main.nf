//
// Run Variant recalibration
//

include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNPS   } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDELS } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNPS   } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDELS } from '../../../modules/nf-core/gatk4/applyvqsr/main'

workflow VCF_VQSR_GATK4 {
    take:
        ch_vcfs         // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_fasta        // channel: [mandatory] [ val(meta2), path(fasta) ]
        ch_fai          // channel: [mandatory] [ val(meta3), path(fai) ]
        ch_dict         // channel: [mandatory] [ val(meta4), path(dict) ]
        ch_hapmap       // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_omni_1000G   // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_snps_1000G   // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_dbsnp        // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_indels_1000G // channel: [mandatory] [ path(vcf), path(tbi) ]

    main:

    ch_versions = Channel.empty()

    ch_hapmap.map { vcf, tbi ->
        def label = "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${vcf.name}"
        [ label, vcf, tbi ]
    }.mix(
        ch_omni_1000G.map { vcf, tbi ->
            def label = "--resource:omni,known=false,training=true,truth=false,prior=12.0 ${vcf.name}"
            [ label, vcf, tbi ]
        },
        ch_snps_1000G.map { vcf, tbi ->
            def label = "--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${vcf.name}"
            [ label, vcf, tbi ]
        },
        ch_dbsnp.map { vcf, tbi ->
            def label = "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${vcf.name}"
            [ label, vcf, tbi ]
        }
    ).multiMap { label, vcf, tbi ->
        labels: [ label ]
        resources: [ vcf, tbi ]
    }
    .set { ch_snp_resources }

    ch_indels_1000G.map { vcf, tbi ->
        def label = "--resource:1000G,known=false,truth=true,training=true,prior=10.0 ${vcf.name}"
        [ label, vcf, tbi ]
    }.mix(
        ch_dbsnp.map { vcf, tbi ->
            def label = "--resource:dbsnp,known=true,truth=false,training=false,prior=2.0 ${vcf.name}"
            [ label, vcf, tbi ]
        }
    ).multiMap { label, vcf, tbi ->
        labels: [ label ]
        resources: [ vcf, tbi ]
    }
    .set { ch_indel_resources }

    GATK4_VARIANTRECALIBRATOR_SNPS(
        ch_vcfs,
        ch_snp_resources.resources.collect(),
        [],
        ch_snp_resources.labels.collect().view(),
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] },
        ch_dict.map { it[1] }
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTRECALIBRATOR_SNPS.out.versions.first())

    GATK4_VARIANTRECALIBRATOR_INDELS(
        ch_vcfs,
        ch_indel_resources.resources.collect(),
        [],
        ch_indel_resources.labels.collect().view(),
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] },
        ch_dict.map { it[1] }
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTRECALIBRATOR_INDELS.out.versions.first())

    ch_vcfs
        .join(GATK4_VARIANTRECALIBRATOR_SNPS.out.recal,    failOnDuplicate:true, failOnMismatch:true)
        .join(GATK4_VARIANTRECALIBRATOR_SNPS.out.idx,      failOnDuplicate:true, failOnMismatch:true)
        .join(GATK4_VARIANTRECALIBRATOR_SNPS.out.tranches, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_applyvqsr_snps_input }

    GATK4_APPLYVQSR_SNPS(
        ch_applyvqsr_snps_input,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] },
        ch_dict.map { it[1] }
    )
    ch_versions = ch_versions.mix(GATK4_APPLYVQSR_SNPS.out.versions.first())

    GATK4_APPLYVQSR_SNPS.out.vcf
        .join(GATK4_APPLYVQSR_SNPS.out.tbi,                  failOnDuplicate:true, failOnMismatch:true)
        .join(GATK4_VARIANTRECALIBRATOR_INDELS.out.recal,    failOnDuplicate:true, failOnMismatch:true)
        .join(GATK4_VARIANTRECALIBRATOR_INDELS.out.idx,      failOnDuplicate:true, failOnMismatch:true)
        .join(GATK4_VARIANTRECALIBRATOR_INDELS.out.tranches, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_applyvqsr_indels_input }

    GATK4_APPLYVQSR_INDELS(
        ch_applyvqsr_indels_input,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] },
        ch_dict.map { it[1] }
    )
    ch_versions = ch_versions.mix(GATK4_APPLYVQSR_INDELS.out.versions.first())

    ch_vcfs = GATK4_APPLYVQSR_INDELS.out.vcf
        .join(GATK4_APPLYVQSR_INDELS.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vcfs = ch_vcfs
    versions = ch_versions     // [ path(versions) ]
}
