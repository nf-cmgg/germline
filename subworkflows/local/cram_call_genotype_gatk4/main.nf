//
// Call and genotype variants with GATK4 tooling
//

include { CRAM_CALL_GATK4               } from '../cram_call_gatk4/main'
include { GVCF_JOINT_GENOTYPE_GATK4     } from '../gvcf_joint_genotype_gatk4/main'
include { VCF_FILTER_BCFTOOLS           } from '../vcf_filter_bcftools/main'
include { VCF_EXTRACT_RELATE_SOMALIER   } from '../vcf_extract_relate_somalier/main'
include { VCF_PED_RTGTOOLS              } from '../vcf_ped_rtgtools/main'

workflow CRAM_CALL_GENOTYPE_GATK4 {
    take:
        ch_input             // channel: [mandatory] [ val(meta), path(cram), path(crai), path(bed) ] => sample CRAM files and their indexes with the split bed files
        ch_gvcfs             // channel: [mandatory] [ val(meta), path(gvcf), path(tbi) ] => earlier called GVCFs with their indices
        ch_peds              // channel: [mandatory] [ val(meta), path(ped) ] => bed files created with the sample preparation subworkflow
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dict              // channel: [mandatory] [ val(meta), path(dict) ] => sequence dictionary
        ch_strtablefile      // channel: [optional]  [ path(strtablefile) ] => STR table file
        ch_dbsnp             // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi         // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF
        ch_somalier_sites    // channel: [optional]  [ path(somalier_sites_vcf) ] => The VCF containing the somalier sites

    main:

    ch_versions     = Channel.empty()
    ch_vcfs         = Channel.empty()
    ch_somalier_ped = Channel.empty()
    ch_reports      = Channel.empty()

    CRAM_CALL_GATK4(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_strtablefile,
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(CRAM_CALL_GATK4.out.versions)
    ch_reports  = ch_reports.mix(CRAM_CALL_GATK4.out.reports)

    ch_gvcfs_ready = ch_gvcfs
        .map { meta, gvcf, tbi ->
            def new_meta = meta + [caller:"haplotypecaller"]
            [ new_meta, gvcf, tbi ]
        }
        .mix(CRAM_CALL_GATK4.out.gvcfs)

    if(!params.only_call) {

        GVCF_JOINT_GENOTYPE_GATK4(
            ch_gvcfs_ready,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_versions = ch_versions.mix(GVCF_JOINT_GENOTYPE_GATK4.out.versions)
        ch_reports  = ch_reports.mix(GVCF_JOINT_GENOTYPE_GATK4.out.reports)

    }

    if(!params.only_call && !params.only_merge) {

        if(params.filter) {
            VCF_FILTER_BCFTOOLS(
                GVCF_JOINT_GENOTYPE_GATK4.out.vcfs
            )
            ch_versions = ch_versions.mix(VCF_FILTER_BCFTOOLS.out.versions)

            VCF_FILTER_BCFTOOLS.out.vcfs
                .set { ch_filter_vcfs }
        } else {
            GVCF_JOINT_GENOTYPE_GATK4.out.vcfs
                .set { ch_filter_vcfs }
        }

        //
        // Run relation tests with somalier
        //

        VCF_EXTRACT_RELATE_SOMALIER(
            ch_filter_vcfs
                .filter { meta, vcf, tbi ->
                    // Filter out the families that only have one individual or are single-sample VCFs
                    meta.family_count > 1
                },
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] },
            ch_somalier_sites,
            ch_peds
                .map { meta, ped ->
                    def new_meta = meta + [caller:"haplotypecaller"]
                    [ new_meta, ped ]
                }
                .filter { meta, ped ->
                    // Filter out the families that only have one individual
                    meta.family_count > 1
                }
        )
        ch_versions = ch_versions.mix(VCF_EXTRACT_RELATE_SOMALIER.out.versions)

        VCF_EXTRACT_RELATE_SOMALIER.out.peds
            .set { ch_somalier_ped }

        //
        // Add PED headers to the VCFs
        //

        if(params.add_ped){
            ch_filter_vcfs
                .branch { meta, vcf, tbi=[] ->
                    // Only add ped headers to VCFs with more than one individual
                    single: meta.family_count == 1
                        return [ meta, vcf ]
                    multiple: meta.family_count > 1
                        return [ meta, vcf ]
                }
                .set { ch_ped_header_branch }

            VCF_PED_RTGTOOLS(
                ch_ped_header_branch.multiple,
                ch_somalier_ped
            )
            ch_versions = ch_versions.mix(VCF_PED_RTGTOOLS.out.versions)

            VCF_PED_RTGTOOLS.out.ped_vcfs
                .mix(ch_ped_header_branch.single)
                .dump(tag:'ped_vcfs', pretty:true)
                .set { ch_vcfs }
        } else {
            ch_filter_vcfs
                .map { meta, vcf, tbi=[] ->
                    [ meta, vcf ]
                }
                .set { ch_vcfs }
        }

    }

    emit:
    vcfs = ch_vcfs         // channel: [ val(meta), path(vcf), path(tbi) ]
    peds = ch_somalier_ped // channel: [ val(meta), path(tsv) ]
    
    reports = ch_reports   // channel: [ path(reports) ]
    versions = ch_versions // channel: [ versions.yml ]

}