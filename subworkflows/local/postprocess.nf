//
// GENOTYPE
//

include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS     } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_COMBINEGVCFS as COMBINEGVCFS        } from '../../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { GATK4_REBLOCKGVCF as REBLOCKGVCF          } from '../../modules/nf-core/modules/gatk4/reblockgvcf/main'
include { TABIX_TABIX as TABIX_GVCFS                } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_COMBINED_GVCFS       } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIP as BGZIP_GENOTYPED_VCFS       } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_PED_VCFS  } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { RTGTOOLS_PEDFILTER as PEDFILTER           } from '../../modules/local/rtgtools/pedfilter/main'
include { MERGE_VCF_HEADERS                         } from '../../modules/local/merge_vcf_headers'
include { BCFTOOLS_FILTER as FILTER_SNPS            } from '../../modules/nf-core/modules/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_INDELS          } from '../../modules/nf-core/modules/bcftools/filter/main'
include { BCFTOOLS_MERGE                            } from '../../modules/nf-core/modules/bcftools/merge/main'
include { BCFTOOLS_CONVERT                          } from '../../modules/nf-core/modules/bcftools/convert/main'
include { BCFTOOLS_VIEW                             } from '../../modules/nf-core/modules/bcftools/view/main'

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

    indexed_gvcfs = gvcfs
                    .combine(TABIX_GVCFS.out.tbi, by: 0)
                    .map({ meta, gvcf, tbi ->
                        [ meta, gvcf, tbi, []]
                    })

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

    combine_gvcfs_input = REBLOCKGVCF.out.vcf
                            .map({ meta, gvcf, tbi ->
                                def new_meta = [:]
                                new_meta.id = meta.family
                                new_meta.family = meta.family

                                [ new_meta, gvcf, tbi ]
                            })
                            .groupTuple()

    //
    // Merge/Combine all the GVCFs from each family
    //

    if (use_bcftools_merge){

        BCFTOOLS_MERGE(
            combine_gvcfs_input,
            [],
            fasta,
            fasta_fai
        )

        combined_gvcfs = BCFTOOLS_MERGE.out.merged_variants
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    } else {

        COMBINEGVCFS(
            combine_gvcfs_input,
            fasta,
            fasta_fai,
            dict
        )

        combined_gvcfs = COMBINEGVCFS.out.combined_gvcf
        ch_versions = ch_versions.mix(COMBINEGVCFS.out.versions)

    }

    //
    // Create indexes for the combined GVCFs
    //

    TABIX_COMBINED_GVCFS(
        combined_gvcfs
    )

    ch_versions = ch_versions.mix(TABIX_COMBINED_GVCFS.out.versions)

    indexed_combined_gvcfs = combined_gvcfs
                            .combine(TABIX_COMBINED_GVCFS.out.tbi, by:0)

    if (!skip_genotyping){

        //
        // Genotype the combined GVCFs
        //

        genotype_gvcfs_input = indexed_combined_gvcfs
                            .map({ meta, gvcf, tbi ->
                                [ meta, gvcf, tbi, [], [] ]
                            })

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

        converted_vcfs = BGZIP_GENOTYPED_VCFS.out.output
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

        converted_vcfs = BCFTOOLS_CONVERT.out.vcf
        ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions)
    }

    //
    // Add pedigree information
    //

    PEDFILTER(
        peds
    )

    ch_versions = ch_versions.mix(PEDFILTER.out.versions)

    merge_vcf_headers_input = converted_vcfs
                                .combine(PEDFILTER.out.vcf, by:0)

    MERGE_VCF_HEADERS(
        merge_vcf_headers_input
    )

    BGZIP_TABIX_PED_VCFS(
        MERGE_VCF_HEADERS.out.vcf
    )

    ch_versions = ch_versions.mix(MERGE_VCF_HEADERS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_TABIX_PED_VCFS.out.versions)

    vcfs_without_index = BGZIP_TABIX_PED_VCFS.out.gz_tbi.map({ meta, vcf, tbi -> [ meta, vcf ]})

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

        post_processed_vcfs = FILTER_INDELS.out.vcf
    }
    else {
        post_processed_vcfs = vcfs_without_index
    }

    emit:
    post_processed_vcfs
    versions = ch_versions
}
