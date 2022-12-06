//
// ANNOTATION
//

include { ENSEMBLVEP                          } from '../../modules/nf-core/ensemblvep/main'
include { VCFANNO                             } from '../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIP as BGZIP_ANNOTATED_VCFS } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_ENSEMBLVEP     } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_CONCAT                     } from '../../modules/nf-core/bcftools/concat/main'

workflow ANNOTATION {
    take:
        vcfs                 // channel: [mandatory] [ meta, vcfs ] => The post-processed VCFs
        fasta                // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai            // channel: [mandatory] [ fasta_fai ] => fasta index
        genome               // value:   [mandatory] Which genome was used to align the samples to
        species              // value:   [mandatory] Which species the samples are from
        vep_cache_version    // value:   [mandatory] which version of VEP to use
        vep_cache            // channel: [optional]  [ vep_cache ] => The VEP cache to use
        vep_extra_files      // channel: [optional]  [ file_1, file_2, file_3, ... ] => All files necessary for using the desired plugins
        vcfanno              // boolean: [mandatory] Whether or not annotation using VCFanno should be performed too
        vcfanno_config       // channel: [mandatory if vcfanno == true] [ toml_config_file ] => The TOML config file for VCFanno
        vcfanno_lua          // channel: [optional  if vcfanno == true] [ lua_file ] => A VCFanno Lua file
        vcfanno_resources    // channel: [mandatory if vcfanno == true] [ resource_dir ] => The directory containing the reference files for VCFanno

    main:

    ch_annotated_vcfs   = Channel.empty()
    ch_reports          = Channel.empty()
    ch_versions         = Channel.empty()

    //
    // Define the present chromosomes
    //

    fasta_fai
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
        .groupTuple()
        .map(
            { meta, chroms ->
                [ meta, chroms.join(",") ]
            }
        )
        .dump(tag:"all_regions", pretty: true)
        .tap { all_regions }
        .count()
        .dump(tag:"count_chromosomes", pretty:true)
        .set { count_chromosomes }

    vcfs
        .combine(all_regions)
        .map(
            { meta, vcf, meta2, chroms ->
                meta = meta + [id:"${meta.id}_${meta2.chr}", regions:chroms]
                [ meta, vcf ]
            }
        )
        .dump(tag:'vep_input', pretty:true)
        .set { vep_input }

    //
    // Annotate using Ensembl VEP
    //

    ENSEMBLVEP(
        vep_input,
        genome,
        species,
        vep_cache_version,
        vep_cache,
        fasta,
        vep_extra_files
    )

    ch_reports  = ch_reports.mix(ENSEMBLVEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions)

    TABIX_ENSEMBLVEP(
        ENSEMBLVEP.out.vcf
    )

    ch_versions = ch_versions.mix(TABIX_ENSEMBLVEP.out.versions)

    ENSEMBLVEP.out.vcf
        .join(TABIX_ENSEMBLVEP.out.tbi)
        .combine(count_chromosomes)
        .map(
            { meta, vcf, tbi, count ->
                meta = meta + [id:meta.family]
                meta.remove("regions")
                [ groupKey(meta, count), vcf, tbi ]
            }
        )
        .groupTuple()
        .dump(tag:"ensemblvep_concat_input", pretty:true)
        .set { ensemblvep_concat_input }

    BCFTOOLS_CONCAT(
        ensemblvep_concat_input
    )

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    if (vcfanno) {

        BCFTOOLS_CONCAT.out.vcf
            .map(
                { meta, vcf ->
                    [ meta, vcf, [] ]
                }
            )
            .dump(tag:'vcfanno_input', pretty:true)
            .set { vcfanno_input }

        VCFANNO(
            vcfanno_input,
            vcfanno_config,
            vcfanno_lua,
            vcfanno_resources
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
