//
// ANNOTATION
//

include { ENSEMBLVEP_VEP                      } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                             } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIP as BGZIP_ANNOTATED_VCFS } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_ENSEMBLVEP     } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_CONCAT                     } from '../../../modules/nf-core/bcftools/concat/main'

include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF as VCF_ANNOTATE_ENSEMBLVEP } from '../../../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'

workflow VCF_ANNOTATION {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf) ] => The post-processed VCFs
        ch_fasta                // channel: [mandatory] [ val(meta2), path(fasta) ] => fasta reference
        ch_fai                  // channel: [mandatory] [ val(meta3), path(fai) ] => fasta index
        ch_vep_cache            // channel: [optional]  [ path(vep_cache) ] => The VEP cache to use
        ch_vep_extra_files      // channel: [optional]  [ path(file_1, file_2, file_3, ...) ] => All files necessary for using the desired plugins
        ch_vcfanno_config       // channel: [mandatory if params.vcfanno == true] [ path(toml_config_file) ] => The TOML config file for VCFanno
        ch_vcfanno_lua          // channel: [optional  if params.vcfanno == true] [ path(lua_file) ] => A VCFanno Lua file
        ch_vcfanno_resources    // channel: [mandatory if params.vcfanno == true] [ path(resource_dir) ] => The directory containing the reference files for VCFanno
        genome                  // string:  The genome needed for VEP
        species                 // string:  The species to use for VEP
        vep_cache_version       // integer: The version of the VEP cache to use
        vep_chunk_size          // integer: The size of each chunk to split VCFs into
        vcfanno                 // boolean: Use VCFanno for annotation

    main:

    ch_annotated_vcfs   = Channel.empty()
    ch_reports          = Channel.empty()
    ch_versions         = Channel.empty()

    ch_vcfs
        .branch { meta, vcf, tbi=[] ->
            tbi: tbi
            no_tbi: !tbi
                return [ meta, vcf ]
        }
        .set { ch_tabix_input }

    TABIX_ENSEMBLVEP(
        ch_tabix_input.no_tbi
    )
    ch_versions = ch_versions.mix(TABIX_ENSEMBLVEP.out.versions.first())

    ch_tabix_input.no_tbi
        .join(TABIX_ENSEMBLVEP.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_tabix_input.tbi)
        .map { meta, vcf, tbi ->
            [ meta, vcf, tbi, [] ]
        }
        .set { ch_vep_input }

    //
    // Do the VEP annotation
    //

    VCF_ANNOTATE_ENSEMBLVEP(
        ch_vep_input,
        ch_fasta,
        genome == "hg38" ? "GRCh38" : genome,
        species,
        vep_cache_version,
        ch_vep_cache,
        ch_vep_extra_files,
        [], [],
        ["ensemblvep"],
        vep_chunk_size
    )

    ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    ch_reports  = ch_reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vep_reports)

    //
    // Annotate the VCFs with VCFanno
    //

    if (vcfanno) {

        VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi
            .map { meta, vcf, tbi ->
                [ meta, vcf, tbi, [] ]
            }
            .dump(tag:'vcfanno_input', pretty:true)
            .set { ch_vcfanno_input }

        VCFANNO(
            ch_vcfanno_input,
            ch_vcfanno_config,
            ch_vcfanno_lua,
            ch_vcfanno_resources
        )
        ch_versions = ch_versions.mix(VCFANNO.out.versions.first())

        BGZIP_ANNOTATED_VCFS(
            VCFANNO.out.vcf
        )
        ch_versions = ch_versions.mix(BGZIP_ANNOTATED_VCFS.out.versions.first())

        BGZIP_ANNOTATED_VCFS.out.output.set { ch_annotated_vcfs }
    }
    else {
        VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi
            .map { meta, vcf, tbi -> [ meta, vcf ]}
            .set { ch_annotated_vcfs }
    }

    emit:
    annotated_vcfs  = ch_annotated_vcfs   // [ val(meta), path(vcf) ]
    reports         = ch_reports          // [ path(reports) ]
    versions        = ch_versions         // [ path(versions) ]
}
