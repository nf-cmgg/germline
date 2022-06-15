//
// GENOTYPE
//

include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS  } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_COMBINEGVCFS as COMBINEGVCFS     } from '../../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { GATK4_REBLOCKGVCF as REBLOCKGVCF       } from '../../modules/nf-core/modules/gatk4/reblockgvcf/main'
include { TABIX_TABIX as TABIX_GVCFS             } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_COMBINED_GVCFS    } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIP as BGZIP                   } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { RTGTOOLS_PEDFILTER as PEDFILTER        } from '../../modules/nf-core/modules/rtgtools/pedfilter/main'
include { MERGE_VCF_HEADERS                      } from '../../modules/local/merge_vcf_headers'

workflow GENOTYPE {
    take:
        gvcfs        // channel: [mandatory] gvcfs
        peds         // channel: [mandatory] peds
        fasta        // channel: [mandatory] fasta reference
        fasta_fai    // channel: [mandatory] fasta reference index
        dict         // channel: [mandatory] sequence dictionary

    main:

    genotyped_gvcfs  = Channel.empty()
    ch_versions      = Channel.empty()

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

    //
    // Combine the GVCFs in each family
    //

    combine_gvcfs_input = REBLOCKGVCF.out.vcf.map(
    { meta, gvcf, tbi ->
        def new_meta = [:]
        new_meta.id = meta.family

        [ new_meta, gvcf, tbi ]
    })
    .groupTuple()

    COMBINEGVCFS(
        combine_gvcfs_input,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(COMBINEGVCFS.out.versions)

    //
    // Create indexes for the combined GVCFs
    //

    combined_gvcfs = COMBINEGVCFS.out.combined_gvcf

    TABIX_COMBINED_GVCFS(
        combined_gvcfs
    )

    ch_versions = ch_versions.mix(TABIX_COMBINED_GVCFS.out.versions)

    indexed_combined_gvcfs = combined_gvcfs
                            .combine(TABIX_COMBINED_GVCFS.out.tbi, by:0)
                            
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

    BGZIP(
        GENOTYPE_GVCFS.out.vcf
    )

    ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

    //
    // Add pedigree information
    //

    PEDFILTER(
        peds
    )

    merge_vcf_headers_input = BGZIP.out.output
                                .combine(PEDFILTER.out.vcf, by:0)

    MERGE_VCF_HEADERS(
        merge_vcf_headers_input
    )

    genotyped_vcfs = MERGE_VCF_HEADERS.out.vcf

    emit:
    genotyped_vcfs    
    versions = ch_versions
}