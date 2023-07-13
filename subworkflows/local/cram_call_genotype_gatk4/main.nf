//
// Call and genotype variants with GATK4 tooling
//

include { CRAM_CALL_GATK4           } from '../cram_call_gatk4/main'
include { GVCF_JOINT_GENOTPE_GATK4  } from '../gvcf_joint_genotype_gatk4/main'
include { VCF_FILTER_BCFTOOLS       } from '../vcf_filter_bcftools/main'

workflow CRAM_CALL_GENOTYPE_GATK4 {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their indexes
        ch_gvcfs             // channel: [mandatory] [ val(meta), path(gvcf), path(tbi) ] => earlier called GVCFs with their indices
        ch_beds              // channel: [mandatory] [ val(meta), path(bed) ] => bed files created with the sample preparation subworkflow
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dict              // channel: [mandatory] [ val(meta), path(dict) ] => sequence dictionary
        ch_strtablefile      // channel: [optional]  [ path(strtablefile) ] => STR table file
        ch_dbsnp             // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi         // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF

    main:

    ch_versions  = Channel.empty()
    ch_vcfs      = Channel.empty()

    CRAM_CALL_GATK4(
        ch_crams,
        ch_beds,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_strtablefile,
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(CRAM_CALL_GATK4.out.versions)

    ch_gvcfs = ch_gvcfs.mix(CRAM_CALL_GATK4.out.gvcfs)

    if(!params.only_call) {

        GVCF_JOINT_GENOTPE_GATK4(
            ch_gvcfs,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_versions = ch_versions.mix(GVCF_JOINT_GENOTPE_GATK4.out.versions)

        if(params.filter) {
            VCF_FILTER_BCFTOOLS(
                GVCF_JOINT_GENOTPE_GATK4.out.vcfs
            )
            ch_versions = ch_versions.mix(VCF_FILTER_BCFTOOLS.out.versions)

            VCF_FILTER_BCFTOOLS.out.vcfs
                .set { ch_vcfs }
        } else {
            GVCF_JOINT_GENOTPE_GATK4.out.vcfs
                .set { ch_vcfs }
        }
    }

    emit:
    vcfs = ch_vcfs         // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions // channel: [ versions.yml ]

}