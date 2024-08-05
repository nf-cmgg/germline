//
// Call and genotype variants with GATK4 tooling
//

include { CRAM_CALL_GATK4               } from '../cram_call_gatk4/main'
include { GVCF_JOINT_GENOTYPE_GATK4     } from '../gvcf_joint_genotype_gatk4/main'
include { VCF_VQSR_GATK4                } from '../vcf_vqsr_gatk4/main'
include { VCF_FILTER_BCFTOOLS           } from '../vcf_filter_bcftools/main'

workflow CRAM_CALL_GENOTYPE_GATK4 {
    take:
        ch_input            // channel: [mandatory] [ val(meta), path(cram), path(crai), path(bed) ] => sample CRAM files and their indexes with the split bed files
        ch_gvcfs            // channel: [mandatory] [ val(meta), path(gvcf), path(tbi) ] => earlier called GVCFs with their indices
        ch_fasta            // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai              // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dict             // channel: [mandatory] [ val(meta), path(dict) ] => sequence dictionary
        ch_strtablefile     // channel: [optional]  [ path(strtablefile) ] => STR table file
        ch_dbsnp            // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi        // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF
        ch_hapmap           // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_omni_1000G       // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_snps_1000G       // channel: [mandatory] [ path(vcf), path(tbi) ]
        ch_indels_1000G     // channel: [mandatory] [ path(vcf), path(tbi) ]
        dragstr             // boolean: create a DragSTR model and run haplotypecaller with it
        only_call           // boolean: only run the variant calling
        only_merge          // boolean: run until the family merging
        vqsr                // boolean: run variant recalibration
        filter              // boolean: filter the VCFs
        scatter_count       // integer: the amount of times the VCFs should be scattered

    main:

    ch_versions     = Channel.empty()
    ch_vcfs         = Channel.empty()
    ch_reports      = Channel.empty()

    CRAM_CALL_GATK4(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_strtablefile,
        ch_dbsnp,
        ch_dbsnp_tbi,
        dragstr
    )
    ch_versions = ch_versions.mix(CRAM_CALL_GATK4.out.versions)
    ch_reports  = ch_reports.mix(CRAM_CALL_GATK4.out.reports)

    ch_gvcfs_ready = ch_gvcfs
        .map { meta, gvcf, tbi ->
            def new_meta = meta + [caller:"haplotypecaller"]
            [ new_meta, gvcf, tbi ]
        }
        .mix(CRAM_CALL_GATK4.out.gvcfs)

    if(!only_call) {

        GVCF_JOINT_GENOTYPE_GATK4(
            ch_gvcfs_ready,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi,
            only_merge,
            scatter_count
        )
        ch_versions = ch_versions.mix(GVCF_JOINT_GENOTYPE_GATK4.out.versions)

    }

    if(!only_call && !only_merge) {

        if(vqsr) {
            VCF_VQSR_GATK4(
                GVCF_JOINT_GENOTYPE_GATK4.out.vcfs,
                ch_fasta,
                ch_fai,
                ch_dict,
                ch_hapmap,
                ch_omni_1000G,
                ch_snps_1000G,
                ch_dbsnp.combine(ch_dbsnp_tbi).collect(),
                ch_indels_1000G,
            )
            ch_versions = ch_versions.mix(VCF_VQSR_GATK4.out.versions)

            VCF_VQSR_GATK4.out.vcfs.set { ch_vqsr_vcfs }

        } else {
            GVCF_JOINT_GENOTYPE_GATK4.out.vcfs
                .set { ch_vqsr_vcfs }
        }

        if(filter) {
            VCF_FILTER_BCFTOOLS(
                ch_vqsr_vcfs,
                true
            )
            ch_versions = ch_versions.mix(VCF_FILTER_BCFTOOLS.out.versions)

            VCF_FILTER_BCFTOOLS.out.vcfs
                .set { ch_vcfs }
        } else {
            ch_vqsr_vcfs
                .set { ch_vcfs }
        }

    }

    emit:
    vcfs = ch_vcfs         // channel: [ val(meta), path(vcf), path(tbi) ]

    reports = ch_reports   // channel: [ path(reports) ]
    versions = ch_versions // channel: [ versions.yml ]

}
