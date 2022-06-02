//
// GENOTYPE
//

include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS  } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_COMBINEGVCFS as COMBINEGVCFS     } from '../../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { TABIX_TABIX as TABIX_GVCFS             } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_COMBINED_GVCFS    } from '../../modules/nf-core/modules/tabix/tabix/main'

workflow GENOTYPE {
    take:
        gvcfs   // channel: [mandatory] gvcfs

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

    ch_versions = ch_versions.mix(TABIX_GVCFS.out.versions)

    //
    // Combine the GVCFs in each family
    //

    combine_gvcfs_input = indexed_gvcfs.map(
    { meta, gvcf, tbi ->
        def new_meta = [:]
        new_meta.id = meta.family

        [ new_meta, gvcf, tbi ]
    })
    .groupTuple()

    COMBINEGVCFS(
        combine_gvcfs_input,
        params.fasta,
        params.fasta_fai,
        params.dict
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
        params.fasta,
        params.fasta_fai,
        params.dict,
        [],
        []
    )

    ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

    genotyped_gvcfs = GENOTYPE_GVCFS.out.vcf

    emit:
    genotyped_gvcfs    
    versions = ch_versions
}