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
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf) ] => The post-processed VCFs
        ch_fasta                // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai                  // channel: [mandatory] [ path(fai) ] => fasta index
        ch_vep_cache            // channel: [optional]  [ path(vep_cache) ] => The VEP cache to use
        ch_vep_extra_files      // channel: [optional]  [ path(file_1, file_2, file_3, ...) ] => All files necessary for using the desired plugins
        ch_vcfanno_config       // channel: [mandatory if params.vcfanno == true] [ path(toml_config_file) ] => The TOML config file for VCFanno
        ch_vcfanno_lua          // channel: [optional  if params.vcfanno == true] [ path(lua_file) ] => A VCFanno Lua file
        ch_vcfanno_resources    // channel: [mandatory if params.vcfanno == true] [ path(resource_dir) ] => The directory containing the reference files for VCFanno

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
                // If the chromosome is a full chromosome, add it as the chr meta
                meta = [chr:chrom]
            }
            else if (chrom ==~ /^.*_alt.*$/) {
                // If the chromosome is an alt contig, add alt as the chr meta
                meta = [chr:"alt"]
            }
            else {
                // If the chromosome isn't a full chromosome or an alt contig, add other as the chr meta
                meta = [chr:"other"]
            }
            [ meta, chrom ]
        }
        .filter { meta, chrom ->
            // Filter out the chromosomes tagged as other
            meta.chr != "other"
        }
        .groupTuple() // No size needed here since this always originates from one file
        .map { meta, chroms ->
            // Join all alt contigs by commas (these will together be used as intervals)
            [ meta, chroms.join(",") ]
        }
        .dump(tag:"all_regions", pretty: true)
        .tap { ch_all_regions }
        .count() // Count how much chromosomes are found (alt contigs are counted as one for all contigs)
        .dump(tag:"count_chromosomes", pretty:true)
        .set { ch_count_chromosomes }

    ch_vcfs
        .combine(ch_all_regions)
        .map { meta, vcf, meta2, chroms ->
            // Remove the chromosomes meta and add the found chromosomes to the VCF meta
            [ meta + [regions:chroms], vcf, [] ]
        }
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
        ch_fasta.map { [[], it] },
        ch_vep_extra_files
    )
    ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())

    ENSEMBLVEP_VEP.out.vcf
        .dump(tag:'vep_output', pretty:true)
        .set { ch_vep_vcfs }

    //
    // Create the VCF index for the annotated VCFs and concatenate all VCFs from the same family back together
    //

    TABIX_ENSEMBLVEP(
        ch_vep_vcfs
    )
    ch_versions = ch_versions.mix(TABIX_ENSEMBLVEP.out.versions.first())

    ch_vep_vcfs
        .join(TABIX_ENSEMBLVEP.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .combine(ch_count_chromosomes)
        .map { meta, vcf, tbi, count ->
            // Remove the chromosome regions from the meta and specify the group size as the amount of chromosomes found earlier
            new_meta = meta - meta.subMap("regions")
            [ groupKey(new_meta, count), vcf, tbi ]
        }
        .groupTuple()
        .dump(tag:"ensemblvep_concat_input", pretty:true)
        .set { ch_ensemblvep_concat_input }

    BCFTOOLS_CONCAT(
        ch_ensemblvep_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    //
    // Annotate the VCFs with VCFanno
    //

    if (params.vcfanno) {

        BCFTOOLS_CONCAT.out.vcf
            .map { meta, vcf ->
                [ meta, vcf, [], [] ]
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
        BCFTOOLS_CONCAT.out.vcf.set { ch_annotated_vcfs }
    }


    emit:
    annotated_vcfs  = ch_annotated_vcfs // [ val(meta), path(vcf) ]
    reports         = ch_reports        // [ path(reports) ]
    versions        = ch_versions       // [ path(versions) ]
}
