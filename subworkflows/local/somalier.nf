//
// SOMALIER
//

include { SOMALIER_EXTRACT } from '../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../modules/nf-core/somalier/relate/main'
include { TABIX_TABIX      } from '../../modules/nf-core/tabix/tabix/main'

workflow SOMALIER {
    take:
        vcfs                 // channel: [mandatory] [ meta, vcfs ] => The joint-genotyped VCFs
        fasta                // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai            // channel: [mandatory] [ fai ] => index of the fasta reference
        somalier_sites       // channel: [mandatory] [ somalier_sites_vcf ] => The VCF containing the common sites for Somalier
        peds                 // channel: [mandatory] [ meta, ped ] => The peds with their meta fields
    main:

    ch_versions         = Channel.empty()

    TABIX_TABIX(
        vcfs
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    vcfs
        .join(TABIX_TABIX.out.tbi)
        .dump(tag:'somalierextract_input', pretty:true)
        .set { somalierextract_input }

    SOMALIER_EXTRACT(
        somalierextract_input,
        fasta,
        fasta_fai,
        somalier_sites
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions)

    SOMALIER_EXTRACT.out.extract
        .join(peds)
        .dump(tag:'somalierrelate_input', pretty:true)
        .set { somalierrelate_input }

    SOMALIER_RELATE(
        somalierrelate_input,
        []
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions)

    SOMALIER_RELATE.out.samples_tsv
        .dump(tag:'generated_peds', pretty:true)
        .set { generated_peds }

    emit:
    generated_peds
    versions        = ch_versions
}
