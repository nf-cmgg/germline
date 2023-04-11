//
// ANNOTATION
//

include { ENSEMBLVEP_VEP                      } from '../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                             } from '../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIP as BGZIP_ANNOTATED_VCFS } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_ENSEMBLVEP     } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_CONCAT                     } from '../../modules/nf-core/bcftools/concat/main'

workflow ANNOTATION {
    take:
        ch_vcfs                 // channel: [mandatory] [ meta, vcfs ] => The post-processed VCFs
        ch_fasta                // channel: [mandatory] [ fasta ] => fasta reference
        ch_fai                  // channel: [mandatory] [ fasta_fai ] => fasta index
        ch_vep_cache            // channel: [optional]  [ vep_cache ] => The VEP cache to use
        ch_vep_extra_files      // channel: [optional]  [ file_1, file_2, file_3, ... ] => All files necessary for using the desired plugins
        ch_vcfanno_config       // channel: [mandatory if params.vcfanno == true] [ toml_config_file ] => The TOML config file for VCFanno
        ch_vcfanno_lua          // channel: [optional  if params.vcfanno == true] [ lua_file ] => A VCFanno Lua file
        ch_vcfanno_resources    // channel: [mandatory if params.vcfanno == true] [ resource_dir ] => The directory containing the reference files for VCFanno

    main:

    ch_annotated_vcfs   = Channel.empty()
    ch_reports          = Channel.empty()
    ch_versions         = Channel.empty()

    //
    // Define the present chromosomes
    //

    ch_fai
        .splitText() { region ->
                chrom = (region[0].tokenize("\t")[0])
                if (chrom ==~ /^(chr)?[0-9XY]{1,2}$/) {
                    meta = [chr:chrom]
                }
                else if (chrom ==~ /^.*_alt.*$/) {
                    meta = [chr:"alt"]
                }
                else {
                    meta = [chr:"other"]
                }

                [ meta, chrom ]
            }
        .filter(
            { meta, chrom ->
                meta.chr != "other"
            }
        )
        .groupTuple() // No size needed here since this always originates from one file
        .map(
            { meta, chroms ->
                [ meta, chroms.join(",") ]
            }
        )
        .dump(tag:"all_regions", pretty: true)
        .tap { ch_all_regions }
        .count()
        .dump(tag:"count_chromosomes", pretty:true)
        .set { ch_count_chromosomes }

    ch_vcfs
        .combine(ch_all_regions)
        .map(
            { meta, vcf, meta2, chroms ->
                [ meta + [regions:chroms], vcf ]
            }
        )
        .dump(tag:'vep_input', pretty:true)
        .set { ch_vep_input }

    //
    // Annotate using Ensembl VEP
    //

    ENSEMBLVEP_VEP(
        ch_vep_input,
        params.genome,
        params.species,
        params.vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files
    )

    ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    ENSEMBLVEP_VEP.out.vcf
        .dump(tag:'vep_output', pretty:true)
        .set { ch_vep_vcfs }

    TABIX_ENSEMBLVEP(
        ch_vep_vcfs
    )

    ch_versions = ch_versions.mix(TABIX_ENSEMBLVEP.out.versions)

    ch_vep_vcfs
        .join(TABIX_ENSEMBLVEP.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .combine(ch_count_chromosomes)
        .map(
            { meta, vcf, tbi, count ->
                new_meta = meta - meta.subMap("regions")
                [ groupKey(new_meta, count), vcf, tbi ]
            }
        )
        .groupTuple()
        .dump(tag:"ensemblvep_concat_input", pretty:true)
        .set { ch_ensemblvep_concat_input }

    BCFTOOLS_CONCAT(
        ch_ensemblvep_concat_input
    )

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    if (params.vcfanno) {

        BCFTOOLS_CONCAT.out.vcf
            .map { it + [[]]}
            .dump(tag:'vcfanno_input', pretty:true)
            .set { ch_vcfanno_input }

        VCFANNO(
            ch_vcfanno_input,
            ch_vcfanno_config,
            ch_vcfanno_lua,
            ch_vcfanno_resources
        )

        BGZIP_ANNOTATED_VCFS(
            VCFANNO.out.vcf
        )

        ch_annotated_vcfs = BGZIP_ANNOTATED_VCFS.out.output
        ch_versions       = ch_versions.mix(VCFANNO.out.versions)
        ch_versions       = ch_versions.mix(BGZIP_ANNOTATED_VCFS.out.versions)
    }
    else {
        ch_annotated_vcfs = BCFTOOLS_CONCAT.out.vcf
    }


    emit:
    annotated_vcfs  = ch_annotated_vcfs
    reports         = ch_reports
    versions        = ch_versions
}
