include { SAMTOOLS_CONVERT      } from '../../../modules/nf-core/samtools/convert/main'
include { VARDICTJAVA           } from '../../../modules/nf-core/vardictjava/main'
include { TABIX_BGZIPTABIX      } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_STATS        } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_CONCAT_BCFTOOLS           } from '../vcf_concat_bcftools/main'

workflow CRAM_CALL_VARDICTJAVA {
    take:
        ch_input             // channel: [mandatory] [ val(meta), path(cram), path(crai), path(bed) ] => sample CRAM files and their indexes
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index

    main:
        ch_versions = Channel.empty()
        ch_reports  = Channel.empty()

        ch_input
            .map { meta, cram, crai, bed ->
                def new_meta = meta + [id:meta.sample]
                [ new_meta, cram, crai, bed ]
            }
            .tap { ch_original }
            .groupTuple()
            .map { meta, cram, crai, beds ->
                [ meta, cram.unique()[0], crai.unique()[0] ]
            }
            .set { ch_crams }

        ch_crams
            .branch { meta, cram, crai ->
                bam: cram.extension == "bam"
                cram: cram.extension == "cram"
            }
            .set { ch_cram_bam }

        SAMTOOLS_CONVERT(
            ch_cram_bam.cram,
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] }
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

        ch_cram_bam.bam
            .mix(SAMTOOLS_CONVERT.out.alignment_index)
            .combine(ch_original, by:0)
            .map { meta, bam, bai, cram, crai, bed ->
                def new_meta = meta + [id:bed.baseName]
                [ new_meta, bam, bai, bed ]
            }
            .set { ch_vardict_input }

        VARDICTJAVA(
            ch_vardict_input,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

        TABIX_BGZIPTABIX(
            VARDICTJAVA.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())
        
        VCF_CONCAT_BCFTOOLS(
            TABIX_BGZIPTABIX.out.gz_tbi
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

        BCFTOOLS_STATS(
            VCF_CONCAT_BCFTOOLS.out.vcfs,
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
        ch_reports = ch_reports.mix(BCFTOOLS_STATS.out.stats.collect { it[1] })

    emit:
    vcfs = VCF_CONCAT_BCFTOOLS.out.vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    reports = ch_reports                // channel: [ path(reports) ]    
    versions = ch_versions              // channel: [ path(versions.yml) ]

}