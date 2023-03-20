//
// GENOTYPE
//

include { MERGE_BEDS                                 } from '../../modules/local/merge_beds'
include { SPLIT_BEDS                                 } from '../../modules/local/split_beds/main'

include { GATK4_GENOMICSDBIMPORT as GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS      } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BEDTOOLS_MAKEWINDOWS                       } from '../../modules/nf-core/bedtools/makewindows/main'
include { GOLEFT_INDEXSPLIT                          } from '../../modules/nf-core/goleft/indexsplit/main'

include { VCF_GATHER_BCFTOOLS                        } from '../../subworkflows/nf-core/vcf_gather_bcftools/main'

workflow JOINT_GENOTYPING {
    take:
        gvcfs               // channel: [mandatory] [ meta, gvcf, tbi ] => The fresh GVCFs called with HaplotypeCaller
        crams               // channel: [mandatory] [ meta, cram, crai ] => The CRAM files of the individuals
        peds                // channel: [mandatory] [ meta, peds ] => The pedigree files for the samples
        fasta               // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai           // channel: [mandatory] [ fasta_fai ] => fasta reference index
        dict                // channel: [mandatory] [ dict ] => sequence dictionary

    main:

    ch_versions         = Channel.empty()

    crams
        .map { meta, cram, crai ->
                new_meta = [
                    family:         meta.family,
                    id:             meta.family,
                    family_count:   meta.family_count
                ]
                [ groupKey(new_meta, meta.family_count.toInteger()), crai ]
            }
        .groupTuple()
        .dump(tag:'genotype_goleft_input', pretty: true)
        .set { goleft_input }

    GOLEFT_INDEXSPLIT(
        goleft_input,
        fasta_fai.map { [ [], it ]},
        params.scatter_count
    )
    ch_versions = ch_versions.mix(GOLEFT_INDEXSPLIT.out.versions)

    SPLIT_BEDS(
        GOLEFT_INDEXSPLIT.out.bed.map { meta, bed ->
            if(workflow.stubRun){
                (1..params.scatter_count).each {
                    start = 100*it
                    bed << "chr22\t${start}\t${start+50}\t0.5\t1\n"
                }
            }
            [ meta, bed, meta.family_count ]
        }
    )
    ch_versions = ch_versions.mix(SPLIT_BEDS.out.versions.first())

    gvcfs
        .map(
            { meta, gvcf, tbi ->
                new_meta = [
                    family:         meta.family,
                    id:             meta.family,
                    family_count:   meta.family_count
                ]
                [ groupKey(new_meta, meta.family_count.toInteger()), gvcf, tbi ]
            }
        )
        .groupTuple()
        .join(SPLIT_BEDS.out.beds, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, gvcfs, tbis, beds ->
            new_meta = meta + [region_count: beds instanceof ArrayList ? beds.size() : 1]
            [ new_meta, gvcfs, tbis, beds ]
        }
        .transpose(by:3)
        .map { meta, gvcfs, tbis, bed ->
            new_meta = meta + [id:bed.baseName]
            [ new_meta, gvcfs, tbis, bed, [], [] ]
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
        .join(GENOTYPE_GVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi ->
            new_meta = meta - meta.subMap("region")
            [ new_meta, vcf, tbi ]
        }
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

    VCF_GATHER_BCFTOOLS.out.vcf
        .map { meta, vcf ->
            new_meta = meta - meta.subMap("region_count")
            [ new_meta, vcf ]
        }
        .set { genotyped_vcfs }

    emit:
    genotyped_vcfs              // channel: [meta, vcf] => The output channel containing the post processed VCF
    versions = ch_versions
}
