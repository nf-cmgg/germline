//
// ANNOTATION
//

include { ENSEMBLVEP } from '../../modules/nf-core/modules/ensemblvep/main'

workflow ANNOTATION {
    take:
        vcfs              // channel: [mandatory] [meta, vcfs]
        fasta             // channel: [mandatory] fasta
        genome            // channel: [mandatory] genome
        species           // channel: [mandatory] species
        vep_cache_version // channel: [mandatory] vep_version
        vep_merged_cache  // channel: [optional] vep_merged_cache
        vep_extra_files   // channel: [optional] vep_extra_files

    main:

    ch_annotated_vcfs   = Channel.empty()
    ch_reports          = Channel.empty()
    ch_versions         = Channel.empty()

    //
    // Annotate using Ensembl VEP
    //

    ENSEMBLVEP(
        vcfs,
        genome,
        species,
        vep_cache_version,
        vep_merged_cache,
        fasta,
        vep_extra_files
    )

    ch_annotated_gvcfs  = ENSEMBLVEP.out.vcf
    ch_reports          = ch_reports.mix(ENSEMBLVEP.out.report)
    ch_versions         = ch_versions.mix(ENSEMBLVEP.out.versions)


    emit:
    annotated_vcfs  = ch_annotated_vcfs
    reports         = ch_reports 
    versions        = ch_versions
}