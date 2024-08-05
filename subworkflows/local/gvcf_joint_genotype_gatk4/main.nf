//
// GENOTYPE
//

include { MERGE_BEDS             } from '../../../modules/local/merge_beds'

include { GAWK                   } from '../../../modules/nf-core/gawk/main'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BEDTOOLS_SPLIT         } from '../../../modules/nf-core/bedtools/split/main'
include { BCFTOOLS_QUERY         } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_STATS         } from '../../../modules/nf-core/bcftools/stats/main'

include { INPUT_SPLIT_BEDTOOLS   } from '../input_split_bedtools/main'
include { VCF_CONCAT_BCFTOOLS    } from '../vcf_concat_bcftools/main'

workflow GVCF_JOINT_GENOTYPE_GATK4 {
    take:
        ch_gvcfs        // channel: [mandatory] [ val(meta), path(gvcf), path(tbi) ] => The GVCFs called with HaplotypeCaller
        ch_fasta        // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai          // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_dict         // channel: [mandatory] [ path(dict) ] => sequence dictionary
        ch_dbsnp        // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi    // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF
        only_merge      // boolean: Only run until the merging of the VCFs
        scatter_count   // integer: The amount of times each file should be scattered

    main:

    ch_versions = Channel.empty()
    ch_vcfs     = Channel.empty()

    //
    // Get a BED file containing all contigs
    //

    GAWK(
        ch_fai,
        []
    )
    ch_versions = ch_versions.mix(GAWK.out.versions)

    //
    // Create GenomicDBs for each family for each BED file
    //

    ch_gvcfs
        .map { meta, gvcf, tbi ->
            // Create the family meta
            def new_meta = [
                family:         meta.family,
                id:             meta.family,
                family_count:   meta.family_count,
                caller:         meta.caller
            ]
            [ groupKey(new_meta, meta.family_count.toInteger()), gvcf, tbi, meta.sample ]
        }
        .groupTuple()
        .map { meta, gvcf, tbi, samples ->
            def new_meta = meta + [samples: "${samples.join(',')}"] // Having a comma-separated string ensures that joins don't fail
            [ new_meta, gvcf, tbi ]
        }
        .combine(GAWK.out.output.map { it[1] })
        .map { meta, gvcfs, tbis, bed ->
            [ meta, gvcfs, tbis, bed, [], [] ]
        }
        .set { ch_genomicsdbimport_input }

    GATK4_GENOMICSDBIMPORT(
        ch_genomicsdbimport_input,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    if(!only_merge) {

        BCFTOOLS_QUERY(
            ch_gvcfs,
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

        BCFTOOLS_QUERY.out.output
            .map { meta, bed ->
                // Create the family meta
                def new_meta = [
                    family:         meta.family,
                    id:             meta.family,
                    family_count:   meta.family_count,
                    caller:         meta.caller
                ]
                [ groupKey(new_meta, meta.family_count.toInteger()), bed, meta.sample ]
            }
            .groupTuple()
            .map { meta, bed, samples ->
                def new_meta = meta + [samples: "${samples.join(',')}"] // Having a comma-separated string ensures that joins don't fail
                [ new_meta, bed ]
            }
            .dump(tag:'merge_beds_input', pretty: true)
            .set { ch_merge_beds_input }

        MERGE_BEDS(
            ch_merge_beds_input,
            ch_fai
        )
        ch_versions = ch_versions.mix(MERGE_BEDS.out.versions.first())

        //
        // Split BED file into multiple BEDs specified by --scatter_count
        //

        INPUT_SPLIT_BEDTOOLS(
            MERGE_BEDS.out.bed.map { meta, bed ->
                // Multiply the scatter count by the family size to better scatter big families
                [meta, bed, (params.scatter_count * meta.family_count)]
            },
            GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, genomicsdb -> [ meta, genomicsdb, [] ]}.view()
        )
        ch_versions = ch_versions.mix(INPUT_SPLIT_BEDTOOLS.out.versions)

        INPUT_SPLIT_BEDTOOLS.out.split.view()
            .map { meta, genomicsdb, extra, bed ->
                [ meta, genomicsdb, [], bed, [] ]
            }
            .set { ch_genotypegvcfs_input }

        //
        // Genotype the genomicsDBs
        //

        GATK4_GENOTYPEGVCFS(
            ch_genotypegvcfs_input,
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] },
            ch_dict.map { it[1] },
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

        GATK4_GENOTYPEGVCFS.out.vcf
            .join(GATK4_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .set { ch_gather_inputs }

        //
        // Combine the genotyped VCFs from each family back together
        //

        VCF_CONCAT_BCFTOOLS(
            ch_gather_inputs,
            true
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

        VCF_CONCAT_BCFTOOLS.out.vcfs
            .set { ch_vcfs }

    }

    emit:
    vcfs = ch_vcfs         // [ val(meta), path(vcf) ]
    versions = ch_versions // [ path(versions) ]

}
