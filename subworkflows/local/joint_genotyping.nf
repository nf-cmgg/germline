//
// GENOTYPE
//

include { MERGE_BEDS                                 } from '../../modules/local/merge_beds'

include { GATK4_GENOMICSDBIMPORT as GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS      } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { TABIX_BGZIP as BGZIP_GENOTYPED_VCFS        } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_PED_VCFS              } from '../../modules/nf-core/tabix/bgzip/main'
include { BCFTOOLS_CONCAT                            } from '../../modules/nf-core/bcftools/concat/main'

include { BED_SCATTER_GROOVY                         } from '../../subworkflows/local/bed_scatter_groovy'

include { VCF_GATHER_BCFTOOLS                        } from '../../subworkflows/nf-core/vcf_gather_bcftools/main'

workflow JOINT_GENOTYPING {
    take:
        gvcfs               // channel: [mandatory] [ meta, gvcf, tbi ] => The fresh GVCFs called with HaplotypeCaller
        beds                // channel: [mandatory] [ meta, bed ] => The BED files of the individuals
        peds                // channel: [mandatory] [ meta, peds ] => The pedigree files for the samples
        fasta               // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai           // channel: [mandatory] [ fasta_fai ] => fasta reference index
        dict                // channel: [mandatory] [ dict ] => sequence dictionary

    main:

    ch_versions         = Channel.empty()

    beds
        .map(
            { meta, bed ->
                new_meta = [:]
                new_meta.family = meta.family
                new_meta.id = meta.family ?: meta.sample
                new_meta.family_count = meta.family_count
                [ groupKey(new_meta, meta.family_count.toInteger()), bed ]
            }
        )
        .groupTuple()
        .dump(tag:'bedmerge_input', pretty: true)
        .set { bedmerge_input }

    MERGE_BEDS(
        bedmerge_input
    )

    ch_versions = ch_versions.mix(MERGE_BEDS.out.versions)

    MERGE_BEDS.out.bed
        .dump(tag:'merged_beds_genotyping', pretty:true)
        .set { beds_to_scatter }

    //
    // Split the merged BED file
    //

    BED_SCATTER_GROOVY(
        beds_to_scatter,
        params.scatter_size
    )

    ch_versions = ch_versions.mix(BED_SCATTER_GROOVY.out.versions)

    BED_SCATTER_GROOVY.out.scattered
        .tap { genomicsdb_beds }
        .map { meta, bed, bed_count ->
            new_meta = meta.findAll{true}[0] + [id:bed.baseName]
            [ new_meta, bed, bed_count ]
        }
        .tap { genotypegvcfs_beds }
        .dump(tag:'scattered_beds_genotyping', pretty:true)
        .set { gather_beds }

    gvcfs
        .map(
            { meta, gvcf, tbi ->
                new_meta = [
                    family: meta.family,
                    id: meta.family ?: meta.sample,
                    family_count: meta.family_count
                ]
                [ groupKey(new_meta, meta.family_count.toInteger()), gvcf, tbi ]
            }
        )
        .groupTuple()
        .combine(genomicsdb_beds, by:0)
        .map(
            { meta, gvcfs, tbis, bed, bed_count ->
                new_meta = meta.findAll{true}[0] + [id:bed.baseName]
                [ new_meta, gvcfs, tbis, bed, [], [] ]
            }
        )
        .set { genomicsdbimport_input }

    genomicsdbimport_input.dump(tag:'genomicsdbimport_input', pretty:true)

    //
    // Merge/Combine all the GVCFs from each family
    //

    GENOMICSDBIMPORT(
        genomicsdbimport_input,
        false,
        false,
        false
    )

    ch_versions = ch_versions.mix(GENOMICSDBIMPORT.out.versions)

    GENOMICSDBIMPORT.out.genomicsdb
        .join(genotypegvcfs_beds, failOnDuplicate: true, failOnMismatch: true)
        .map(
            { meta, genomic_db, bed, bed_count ->
                [ meta, genomic_db, [], bed, [] ]
            }
        )
        .dump(tag:'genotypegvcfs_input', pretty:true)
        .set { genotypegvcfs_input }

    ch_versions = ch_versions.mix(GENOMICSDBIMPORT.out.versions)

    //
    // Genotype the combined GVCFs
    //

    GENOTYPE_GVCFS(
        genotypegvcfs_input,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

    GENOTYPE_GVCFS.out.vcf
        .join(GENOTYPE_GVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .dump(tag:'genotyped_vcfs', pretty:true)
        .set { genotyped_vcfs }

    VCF_GATHER_BCFTOOLS(
        genotyped_vcfs,
        gather_beds,
        "family",
        true
    )

    ch_versions = ch_versions.mix(VCF_GATHER_BCFTOOLS.out.versions)

    emit:
    genotyped_vcfs = VCF_GATHER_BCFTOOLS.out.vcf    // channel: [meta, vcf] => The output channel containing the post processed VCF
    versions = ch_versions
}
