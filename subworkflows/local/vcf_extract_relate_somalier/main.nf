include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'

workflow VCF_EXTRACT_RELATE_SOMALIER {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
        ch_fasta                // channel: [mandatory] [ path(fasta) ]
        ch_fasta_fai            // channel: [mandatory] [ path(fai) ]
        ch_somalier_sites       // channel: [mandatory] [ path(somalier_sites_vcf) ]
        ch_peds                 // channel: [mandatory] [ val(meta), path(ped) ]

    main:

    ch_versions = Channel.empty()

    SOMALIER_EXTRACT(
        ch_vcfs,
        ch_fasta,
        ch_fasta_fai,
        ch_somalier_sites
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions.first())

    SOMALIER_EXTRACT.out.extract
        .join(ch_peds, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, extract, ped ->
            [ meta, extract, ped ]
        }
        .set { ch_somalierrelate_input }

    SOMALIER_RELATE(
        ch_somalierrelate_input,
        []
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions.first())

    emit:
    extract        = SOMALIER_EXTRACT.out.extract       // channel: [ val(meta), path(extract) ]
    html           = SOMALIER_RELATE.out.html           // channel: [ val(meta), path(html) ]
    pairs_tsv      = SOMALIER_RELATE.out.pairs_tsv      // channel: [ val(meta), path(tsv) ]
    samples_tsv    = SOMALIER_RELATE.out.samples_tsv    // channel: [ val(meta), path(tsv) ]
    peds           = SOMALIER_RELATE.out.ped            // channel: [ val(meta), path(tsv) ]
    versions       = ch_versions                        // channel: [ path(versions.yml) ]
}
