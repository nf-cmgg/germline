//
// Call and genotype variants with GATK4 tooling
//

include { CRAM_CALL_GATK4               } from '../cram_call_gatk4/main'
include { GVCF_JOINT_GENOTYPE_GATK4     } from '../gvcf_joint_genotype_gatk4/main'

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
        dragstr             // boolean: create a DragSTR model and run haplotypecaller with it
        only_call           // boolean: only run the variant calling
        only_merge          // boolean: run until the family merging
        filter              // boolean: filter the VCFs
        scatter_count       // integer: the amount of times the VCFs should be scattered

    main:

    def ch_versions     = Channel.empty()
    def ch_vcfs         = Channel.empty()
    def ch_reports      = Channel.empty()

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

    def ch_gvcfs_ready = ch_gvcfs
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

    emit:
    vcfs = GVCF_JOINT_GENOTYPE_GATK4.out.vcfs   // channel: [ val(meta), path(vcf), path(tbi) ]

    reports = ch_reports                        // channel: [ path(reports) ]
    versions = ch_versions                      // channel: [ versions.yml ]

}
