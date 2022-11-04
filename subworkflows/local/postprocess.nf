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
include { TABIX_BGZIP as BGZIP_GENOTYPED_VCFS       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_PED_VCFS  } from '../../modules/nf-core/tabix/bgziptabix/main'
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
        .set { indexed_gvcfs }

    ch_versions = ch_versions.mix(TABIX_GVCFS.out.versions)

    //
    // Reblock the single sample GVCF files
    //

    REBLOCKGVCF(
        indexed_gvcfs,
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
                def new_meta = [:] // TODO don't create a new meta here
                new_meta.id = meta.family
                new_meta.family = meta.family

                [ new_meta, gvcf, tbi ]
            }
        )
        .groupTuple() // TODO add groupTuple size (size of family)
        .set { combine_gvcfs_input }

    //
    // Merge/Combine all the GVCFs from each family
    //

    // TODO add better support for families containing only one sample

    if (use_bcftools_merge){

        BCFTOOLS_MERGE(
            combine_gvcfs_input,
            [],
            fasta,
            fasta_fai
        )

        BCFTOOLS_MERGE.out.merged_variants.set { combined_gvcfs }
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    } else {

        COMBINEGVCFS(
            combine_gvcfs_input,
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
            .set { genotype_gvcfs_input }

        GENOTYPE_GVCFS(
            genotype_gvcfs_input,
            fasta,
            fasta_fai,
            dict,
            [],
            []
        )

        ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

        BGZIP_GENOTYPED_VCFS(
            GENOTYPE_GVCFS.out.vcf
        )

        BGZIP_GENOTYPED_VCFS.out.output.set { converted_vcfs }
        ch_versions = ch_versions.mix(BGZIP_GENOTYPED_VCFS.out.versions)

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

        BCFTOOLS_CONVERT.out.vcf.set { converted_vcfs }
        ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions)
    }

    //
    // Add pedigree information
    //

    // TODO fix header not the same error when one ped file is supplied in a multi sample family

    converted_vcfs
        .combine(peds, by:0)
        .branch(
            { meta, vcf, ped ->
                has_ped: ped
                no_ped: ped == []
            }
        )
    .set { ped_input }

    PEDFILTER(
        ped_input.has_ped.map({ meta, vcf, ped -> [ meta, ped ]})
    )

    ch_versions = ch_versions.mix(PEDFILTER.out.versions)

    converted_vcfs
        .join(PEDFILTER.out.vcf)
        .set { merge_vcf_headers_input }

    MERGE_VCF_HEADERS(
        merge_vcf_headers_input
    )

    BGZIP_TABIX_PED_VCFS(
        MERGE_VCF_HEADERS.out.vcf.mix(ped_input.no_ped.map({ meta, vcf, ped -> [ meta, vcf ]}))
    )

    ch_versions = ch_versions.mix(MERGE_VCF_HEADERS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_TABIX_PED_VCFS.out.versions)

    BGZIP_TABIX_PED_VCFS.out.gz_tbi
        .map({ meta, vcf, tbi -> [ meta, vcf ]})
        .set { vcfs_without_index }

    //
    // Filter the variants
    //

    if (output_mode == "seqplorer") {
        FILTER_SNPS(
            vcfs_without_index
        )

        FILTER_INDELS(
            FILTER_SNPS.out.vcf
        )

        ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)
        ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

        FILTER_INDELS.out.vcf.set { post_processed_vcfs }
    }
    else {
        vcfs_without_index.set { post_processed_vcfs }
    }

    emit:
    post_processed_vcfs
    versions = ch_versions
}
