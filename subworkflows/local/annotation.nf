//
// ANNOTATION
//

include { ENSEMBLVEP                          } from '../../modules/nf-core/ensemblvep/main'
include { VCFANNO                             } from '../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIP as BGZIP_ANNOTATED_VCFS } from '../../modules/nf-core/tabix/bgzip/main'

workflow ANNOTATION {
    take:
        vcfs                 // channel: [mandatory] [ meta, vcfs ] => The post-processed VCFs
        fasta                // channel: [mandatory] [ fasta ] => fasta reference
        genome               // value:   [mandatory] Which genome was used to align the samples to
        species              // value:   [mandatory] Which species the samples are from
        vep_cache_version    // value:   [mandatory] which version of VEP to use
        vep_merged_cache     // channel: [optional]  [ vep_merged_cache ] => The VEP cache to use
        vep_extra_files      // channel: [optional]  [ file_1, file_2, file_3, ... ] => All files necessary for using the desired plugins
        vcfanno              // boolean: [mandatory] Whether or not annotation using VCFanno should be performed too
        vcfanno_toml         // channel: [mandatory if vcfanno == true] [ toml_config_file ] => The TOML config file for VCFanno
        vcfanno_resources    // channel: [mandatory if vcfanno == true] [ resource_dir ] => The directory containing the reference files for VCFanno

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

    if (vcfanno) {
        VCFANNO(
            ENSEMBLVEP.out.vcf.map({ meta, vcf -> [ meta, vcf, [] ] }),
            vcfanno_toml,
            vcfanno_resources
        )

        ch_annotated_vcfs = VCFANNO.out.vcf
        ch_versions       = ch_versions.mix(VCFANNO.out.versions)
    }
    else {
        ch_annotated_vcfs = ENSEMBLVEP.out.vcf
    }

    BGZIP_ANNOTATED_VCFS(
        ch_annotated_vcfs
    )

    ch_versions = ch_versions.mix(BGZIP_ANNOTATED_VCFS.out.versions)

    emit:
    annotated_vcfs  = BGZIP_ANNOTATED_VCFS.out.output
    reports         = ch_reports
    versions        = ch_versions
}
