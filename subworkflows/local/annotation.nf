//
// ANNOTATION
//

include { ENSEMBLVEP } from '../../modules/nf-core/modules/ensemblvep/main'

workflow ANNOTATION {
    take:
        gvcfs   // channel: [mandatory] gvcfs

    main:

    annotated_gvcfs  = Channel.empty()
    ch_versions      = Channel.empty()

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

    annotated_gvcfs = ENSEMBLVEP.out.vcf

    emit:
    annotated_gvcfs    
    versions = ch_versions
}