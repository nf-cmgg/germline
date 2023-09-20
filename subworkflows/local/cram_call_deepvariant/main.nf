include { DEEPVARIANT      } from '../../../modules/nf-core/deepvariant/main'

workflow CRAM_CALL_DEEPVARIANT {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their indexes
        ch_bed               // channel: [mandatory] [ val(meta), path(bed) ] => The BED file to use
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index

    main:
        ch_versions = Channel.empty()

        ch_crams
            .join(ch_bed, failOnDuplicate:true, failOnMismatch:true)
            .map { meta, cram, crai, bed ->
                def new_meta = meta + [caller:"deepvariant"]
                [ new_meta, cram, crai, bed ]
            }
            .set { ch_deepvariant_input }

        DEEPVARIANT(
            ch_deepvariant_input,
            ch_fasta,
            ch_fai,
            [[],[]]
        )
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

        DEEPVARIANT.out.vcf
            .join(DEEPVARIANT.out.vcf_tbi, failOnDuplicate:true, failOnMismatch:true)
            .view()
            .set { ch_vcfs }

    emit:
    vcfs = ch_vcfs          // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions  // channel: [ path(versions.yml) ]

}