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
        .join(MERGE_BEDS.out.bed, failOnDuplicate: true, failOnMismatch: true)
        .map(
            { meta, gvcfs, tbis, bed ->
                new_meta = meta + [region_count:bed.readLines().size()]
                [ new_meta, gvcfs, tbis, bed, [], [] ]
            }
        )
        .splitText(elem:3)
        .map { meta, gvcfs, tbis, region ->
            region_split = region[0..-2].split("\t")
            region_string = "${region_split[0]}:${region_split[1] as int + 1}-${region_split[2]}"
            new_meta = meta + [id:"${meta.id}_${region_string.replace(":","_").replace("-":"_")}", region:region_string]
            [ new_meta, gvcfs, tbis, [], [], [] ]
        }
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
        .map(
            { meta, genomic_db ->
                [ meta, genomic_db, [], [], [] ]
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
        .join(GENOTYPE_GVCFS.out.tbi)
        .dump(tag:'genotyped_vcfs', pretty:true)
        .set { genotyped_vcfs }

    VCF_GATHER_BCFTOOLS(
        genotyped_vcfs,
        genotyped_vcfs.map { meta, vcf, tbi ->
            [ meta, [], meta.region_count ]
        },
        "family",
        true
    )

    ch_versions = ch_versions.mix(VCF_GATHER_BCFTOOLS.out.versions)

    emit:
    genotyped_vcfs = VCF_GATHER_BCFTOOLS.out.vcf    // channel: [meta, vcf] => The output channel containing the post processed VCF
    versions = ch_versions
}
