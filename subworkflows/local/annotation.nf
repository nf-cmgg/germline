//
// ANNOTATION
//

include { ENSEMBLVEP } from '../../modules/nf-core/modules/ensemblvep/main'

workflow ANNOTATION {
    take:
        gvcfs   // channel: [mandatory] gvcfs

    main:

    ch_annotated_gvcfs  = Channel.empty()
    ch_reports          = Channel.empty()
    ch_versions         = Channel.empty()

    //
    // Annotate using Ensembl VEP
    //

    ENSEMBLVEP(
        gvcfs,
        params.genome,
        "homo_sapiens",
        "105",
        [],
        []
    )

    ch_annotated_gvcfs  = ENSEMBLVEP.out.vcf
    ch_reports          = ch_reports.mix(ENSEMBLVEP.out.report)
    ch_versions         = ch_versions.mix(ENSEMBLVEP.out.versions)


    emit:
    annotated_gvcfs = ch_annotated_gvcfs
    reports         = ch_reports 
    versions        = ch_versions
}