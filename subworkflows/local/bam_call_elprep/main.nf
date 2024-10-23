//
// Call the variants using Elprep
//

include { ELPREP_FILTER         } from '../../../modules/nf-core/elprep/filter/main'
include { BCFTOOLS_STATS        } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_CONCAT_BCFTOOLS   } from '../vcf_concat_bcftools/main'
include { VCF_DBSNP_VCFANNO     } from '../vcf_dbsnp_vcfanno/main'

workflow BAM_CALL_ELPREP {
    take:
        ch_input            // channel: [mandatory] [ val(meta), path(bam), path(bai), path(bed) ] => sample BAM files and their indexes with the split bed files
        ch_elfasta          // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_elsites          // channel: [optional]  [ val(meta), path(elsites) ]
        ch_dbsnp            // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi        // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF

    main:

    def ch_versions  = Channel.empty()

    ELPREP_FILTER(
        ch_input.map { meta, bam, bai, bed ->
            def new_meta = meta + [caller:'elprep']
            [ new_meta, bam, bai, bed, [], [], [] ]
        },
        [[],[]],
        ch_elfasta,
        ch_elsites,
        true, // haplotypecaller
        false,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(ELPREP_FILTER.out.versions.first())
    
    VCF_CONCAT_BCFTOOLS(
        ELPREP_FILTER.out.gvcf,
        true
    )
    ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

    def ch_annotated = Channel.empty()
    if(!(ch_dbsnp instanceof List)) {
        VCF_DBSNP_VCFANNO(
            VCF_CONCAT_BCFTOOLS.out.vcfs,
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_versions = ch_versions.mix(VCF_DBSNP_VCFANNO.out.versions)
        ch_annotated = TABIX_BGZIP.out.output
    } else {
        ch_annotated = VCF_CONCAT_BCFTOOLS.out.vcfs
    }

    BCFTOOLS_STATS(
        ch_annotated,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    def ch_reports = BCFTOOLS_STATS.out.stats.collect{ _meta, report -> report}

    emit:
    gvcfs = ch_annotated    // channel: [ val(meta), path(vcf), path(tbi) ]
    reports = ch_reports    // channel: [ path(stats) ]
    versions = ch_versions  // channel: [ versions.yml ]

}
