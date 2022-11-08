//
// GENOTYPE
//

include { RTGTOOLS_PEDFILTER as PEDFILTER           } from '../../modules/local/rtgtools/pedfilter/main'
include { MERGE_VCF_HEADERS                         } from '../../modules/local/merge_vcf_headers'

include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS     } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_COMBINEGVCFS as COMBINEGVCFS        } from '../../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_REBLOCKGVCF as REBLOCKGVCF          } from '../../modules/nf-core/gatk4/reblockgvcf/main'
include { TABIX_TABIX as TABIX_GVCFS                } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_SORTED_GVCFS         } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_COMBINED_GVCFS       } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_POSTPROCESSED_VCFS   } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP as BGZIP_GENOTYPED_VCFS       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_PED_VCFS             } from '../../modules/nf-core/tabix/bgzip/main'
include { BCFTOOLS_FILTER as FILTER_SNPS            } from '../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_INDELS          } from '../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_MERGE                            } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_CONVERT                          } from '../../modules/nf-core/bcftools/convert/main'
include { BCFTOOLS_VIEW                             } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_SORT                             } from '../../modules/nf-core/bcftools/sort/main'

workflow POST_PROCESS {
    take:
        gvcfs               // channel: [mandatory] [ meta, gvcf ] => The fresh GVCFs called with HaplotypeCaller
        peds                // channel: [mandatory] [ meta, peds ] => The pedigree files for the samples
        fasta               // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai           // channel: [mandatory] [ fasta_fai ] => fasta reference index
        dict                // channel: [mandatory] [ dict ] => sequence dictionary
        output_mode         // value:   [mandatory] whether or not to make the output seqplorer- or seqr-compatible
        skip_genotyping     // boolean: [mandatory] whether or not to skip the genotyping
        use_bcftools_merge  // boolean: [mandatory] whether or not to use bcftools merge instead of CombineGVCFs

    main:

    post_processed_vcfs  = Channel.empty()
    ch_versions          = Channel.empty()

    //
    // Create indexes for all the GVCF files
    //

    TABIX_GVCFS(
        gvcfs
    )

    gvcfs
        .join(TABIX_GVCFS.out.tbi)
        .map(
            { meta, gvcf, tbi ->
                [ meta, gvcf, tbi, []]
            }
        )
        .dump(tag:'reblockgvcf_input', pretty:true)
        .set { reblockgvcf_input }

    ch_versions = ch_versions.mix(TABIX_GVCFS.out.versions)

    //
    // Reblock the single sample GVCF files
    //

    REBLOCKGVCF(
        reblockgvcf_input,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    ch_versions = ch_versions.mix(REBLOCKGVCF.out.versions)

    //
    // Sort the VCF files
    //

    BCFTOOLS_SORT(
        REBLOCKGVCF.out.vcf.map({ meta, vcf, tbi -> [ meta, vcf ]})
    )

    TABIX_SORTED_GVCFS(
        BCFTOOLS_SORT.out.vcf
    )

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_SORTED_GVCFS.out.tbi)
        .map(
            { meta, gvcf, tbi ->
                def new_meta = meta.clone()
                new_meta.id = meta.family ?: meta.sample
                new_meta.remove('sample')
                new_meta.family = meta.family

                [ groupKey(new_meta, new_meta.family_count.toInteger()), gvcf, tbi ]
            }
        )
        .groupTuple()
        .branch(
            { meta, vcfs, tbis ->
                multiple: vcfs.size() > 1
                single:   vcfs.size() == 1
            }
        )
        .set { combine_gvcfs_input }

    combine_gvcfs_input.multiple.dump(tag:'combinegvcfs_input_multiple', pretty:true)
    combine_gvcfs_input.single.dump(tag:'combinegvcfs_input_single', pretty:true)

    //
    // Merge/Combine all the GVCFs from each family
    //

    if (use_bcftools_merge){

        BCFTOOLS_MERGE(
            combine_gvcfs_input.multiple,
            [],
            fasta,
            fasta_fai
        )

        BCFTOOLS_MERGE.out.merged_variants.set { combined_gvcfs }
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    } else {

        COMBINEGVCFS(
            combine_gvcfs_input.multiple,
            fasta,
            fasta_fai,
            dict
        )

        COMBINEGVCFS.out.combined_gvcf.set { combined_gvcfs }
        ch_versions = ch_versions.mix(COMBINEGVCFS.out.versions)

    }

    //
    // Create indexes for the combined GVCFs
    //

    TABIX_COMBINED_GVCFS(
        combined_gvcfs
    )

    ch_versions = ch_versions.mix(TABIX_COMBINED_GVCFS.out.versions)

    combined_gvcfs
        .join(TABIX_COMBINED_GVCFS.out.tbi)
        .dump(tag:'combined_gvcfs', pretty:true)
        .mix(combine_gvcfs_input.single)
        .set { indexed_combined_gvcfs }

    if (!skip_genotyping){

        //
        // Genotype the combined GVCFs
        //

        indexed_combined_gvcfs
            .map(
                { meta, gvcf, tbi ->
                    [ meta, gvcf, tbi, [], [] ]
                }
            )
            .dump(tag:'genotypegvcfs_input', pretty:true)
            .set { genotypegvcfs_input }

        GENOTYPE_GVCFS(
            genotypegvcfs_input,
            fasta,
            fasta_fai,
            dict,
            [],
            []
        )

        GENOTYPE_GVCFS.out.vcf.set { converted_vcfs }
        ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

    } else {

        //
        // Remove the ref blocks from the GVCF
        //

        BCFTOOLS_VIEW(
            indexed_combined_gvcfs,
            [],
            [],
            []
        )

        //
        // Convert all the GVCFs to VCF files
        //

        BCFTOOLS_CONVERT(
            BCFTOOLS_VIEW.out.vcf.map({ meta, vcf -> [ meta, vcf, []]}),
            [],
            fasta
        )

        BCFTOOLS_CONVERT.out.vcf_gz.set { converted_vcfs }
        ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions)
    }

    //
    // Add pedigree information
    //

    // TODO fix header not the same error when one ped file is supplied in a multi sample family

    converted_vcfs
        .join(peds)
        .branch(
            { meta, vcf, ped ->
                has_ped: ped
                    return [ meta, vcf, ped ]
                no_ped: !ped
                    return [ meta, vcf ]
            }
        )
        .set { ped_vcfs }

    ped_vcfs.has_ped
        .tap { pedfilter_input }
        .dump(tag:'ped_vcfs_has_ped', pretty:true)
        .set { merge_header_input }
    ped_vcfs.no_ped.dump(tag:'ped_vcfs_no_ped', pretty:true)

    PEDFILTER(
        pedfilter_input.map({ meta, vcf, ped -> [ meta, ped ]})
    )

    ch_versions = ch_versions.mix(PEDFILTER.out.versions)

    BGZIP_GENOTYPED_VCFS(
        merge_header_input
    )

    ch_versions = ch_versions.mix(BGZIP_GENOTYPED_VCFS.out.versions)

    MERGE_VCF_HEADERS(
        BGZIP_GENOTYPED_VCFS.out.output
            .join(PEDFILTER.out.vcf)
    )

    ch_versions = ch_versions.mix(MERGE_VCF_HEADERS.out.versions)

    BGZIP_PED_VCFS(
        MERGE_VCF_HEADERS.out.vcf.mix(ped_vcfs.no_ped.map({ meta, vcf -> [ meta, vcf ]}))
    )

    ch_versions = ch_versions.mix(BGZIP_PED_VCFS.out.versions)

    ped_vcfs.no_ped
        .mix(BGZIP_PED_VCFS.out.output)
        .dump(tag:'filter_input', pretty:true)
        .set { filter_input }

    //
    // Filter the variants
    //

    if (output_mode == "seqplorer") {
        FILTER_SNPS(
            filter_input
        )

        FILTER_INDELS(
            FILTER_SNPS.out.vcf
        )

        ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)
        ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

        FILTER_INDELS.out.vcf.set { post_processed_vcfs }

    }
    else {
        filter_input.set { post_processed_vcfs }
    }

    emit:
    post_processed_vcfs     // channel: [meta, vcf] => The output channel containing the post processed VCF
    versions = ch_versions
}
