//
// ANNOTATION
//

include { ENSEMBLVEP                          } from '../../modules/nf-core/modules/ensemblvep/main'
include { TABIX_BGZIP as BGZIP_ANNOTATED_VCFS } from '../../modules/nf-core/modules/tabix/bgzip/main'

workflow ANNOTATION {
    take:
        vcfs              // channel: [mandatory] [ meta, vcfs ] => The post-processed VCFs
        fasta             // channel: [mandatory] [ fasta ] => fasta reference
        genome            // value:   [mandatory] Which genome was used to align the samples to
        species           // value:   [mandatory] Which species the samples are from
        vep_cache_version // value:   [mandatory] which version of VEP to use
        vep_merged_cache  // channel: [optional]  [ vep_merged_cache ] => The VEP cache to use
        vep_extra_files   // channel: [optional]  [ file_1, file_2, file_3, ... ] => All files necessary for using the desired plugins

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

    ch_reports          = ch_reports.mix(ENSEMBLVEP.out.report)
    ch_versions         = ch_versions.mix(ENSEMBLVEP.out.versions)

    BGZIP_ANNOTATED_VCFS(
        ENSEMBLVEP.out.vcf
    )

    ch_annotated_vcfs   = BGZIP_ANNOTATED_VCFS.out.output
    ch_versions         = ch_versions.mix(BGZIP_ANNOTATED_VCFS.out.versions)

    emit:
    annotated_vcfs  = ch_annotated_vcfs
    reports         = ch_reports 
    versions        = ch_versions
}