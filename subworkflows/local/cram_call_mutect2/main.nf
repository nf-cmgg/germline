//
// Call the variants using HaplotypeCaller
//

include { GATK4_MUTECT2         } from '../../../modules/nf-core/gatk4/mutect2/main'
include { BCFTOOLS_STATS        } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_CONCAT_BCFTOOLS   } from '../vcf_concat_bcftools/main'

workflow CRAM_CALL_MUTECT2 {
    take:
        ch_input                // channel: [mandatory] [ val(meta), path(cram), path(crai), path(bed) ] => sample CRAM files and their indexes with the split bed files
        ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai                  // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dict                 // channel: [mandatory] [ val(meta), path(dict) ] => sequence dictionary
        ch_dbsnp                // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi            // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF
        ch_panel_of_normals     // channel: [optional]  [ path(panel_of_normals) ] => The VCF containing the panel of normals variants
        ch_panel_of_normals_tbi // channel: [optional]  [ path(panel_of_normals_tbi) ] => The index of the panel of normals VCF

    main:

    GATK4_MUTECT2(
        ch_input.map { meta, cram, crai, bed ->
            def new_meta = meta + [caller:'mutect2']
            tuple(new_meta, cram, crai, bed)
        },
        ch_fasta,
        ch_fai.map { meta, fai -> tuple(meta, fai, [])},
        ch_dict,
        [],
        [],
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi
    )

    def ch_calls = GATK4_MUTECT2.out.vcf
        .join(GATK4_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    VCF_CONCAT_BCFTOOLS(
        ch_calls
    )

    BCFTOOLS_STATS(
        VCF_CONCAT_BCFTOOLS.out.vcfs,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )

    emit:
    vcfs    = VCF_CONCAT_BCFTOOLS.out.vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]
    reports = BCFTOOLS_STATS.out.stats      // channel: [ val(meta), path(stats) ]
}
