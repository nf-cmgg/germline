//
// GENOTYPE
//

include { MERGE_BEDS                                 } from '../../modules/local/merge_beds'

include { GATK4_GENOMICSDBIMPORT as GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS      } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BEDTOOLS_SPLIT                             } from '../../modules/nf-core/bedtools/split/main'
include { BCFTOOLS_QUERY                             } from '../../modules/nf-core/bcftools/query/main'

include { VCF_GATHER_BCFTOOLS                        } from '../../subworkflows/nf-core/vcf_gather_bcftools/main'

workflow JOINT_GENOTYPING {
    take:
        ch_gvcfs               // channel: [mandatory] [ val(meta), path(gvcf), path(tbi) ] => The GVCFs called with HaplotypeCaller
        ch_fasta               // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai                 // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_dict                // channel: [mandatory] [ path(dict) ] => sequence dictionary
        ch_dbsnp               // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi           // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF

    main:

    ch_versions = Channel.empty()
    ch_genotyped_vcfs = Channel.empty()


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
            new_meta = [
                family:         meta.family,
                id:             meta.family,
                family_count:   meta.family_count
            ]
            [ groupKey(new_meta, meta.family_count.toInteger()), bed ]
        }
        .groupTuple()
        .dump(tag:'merge_beds_input', pretty: true)
        .set { ch_merge_beds_input }

    MERGE_BEDS(
        ch_merge_beds_input,
        ch_fai
    )
    ch_versions = ch_versions.mix(MERGE_BEDS.out.versions.first())

    //
    // Create GenomicDBs for each family for each BED file
    //

    ch_gvcfs
        .map { meta, gvcf, tbi ->
            // Create the family meta
            new_meta = [
                family:         meta.family,
                id:             meta.family,
                family_count:   meta.family_count
            ]
            [ groupKey(new_meta, meta.family_count.toInteger()), gvcf, tbi ]
        }
        .groupTuple()
        .join(MERGE_BEDS.out.bed, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcfs, tbis, bed -> 
            [ meta, gvcfs, tbis, bed, [], [] ]
        }
        .set { ch_genomicsdbimport_input }

    GENOMICSDBIMPORT(
        ch_genomicsdbimport_input,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(GENOMICSDBIMPORT.out.versions.first())

    if(!params.only_merge) {

        //
        // Split BED file into multiple BEDs specified by --scatter_count
        //

        BEDTOOLS_SPLIT(
            MERGE_BEDS.out.bed.map { meta, bed ->
                // Multiply the scatter count by the family size to better scatter big families
                [meta, bed, (params.scatter_count * meta.family_count) ]
            }
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())


        GENOMICSDBIMPORT.out.genomicsdb
            .join(BEDTOOLS_SPLIT.out.beds, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, genomicsdb, beds ->
                // Determine the amount of BED files per sample
                bed_is_list = beds instanceof ArrayList
                new_meta = meta + [region_count: bed_is_list ? beds.size() : 1]
                [ new_meta, genomicsdb, bed_is_list ? beds : [beds] ]
            }
            .transpose(by:2) // Create one channel entry for each BED file per family
            .map { meta, genomicsdb, bed ->
                // Set the base name of the BED file as the ID (this will look like sample_id.xxxx, where xxxx are numbers)
                new_meta = meta + [id:bed.baseName]
                [ new_meta, genomicsdb, [], bed, [] ]
            }
            .set { ch_genotypegvcfs_input }

        //
        // Genotype the genomicsDBs
        //

        GENOTYPE_GVCFS(
            ch_genotypegvcfs_input,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions.first())

        GENOTYPE_GVCFS.out.vcf
            .join(GENOTYPE_GVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .dump(tag:'gather_inputs_joint_genotyping', pretty:true)
            .set { ch_gather_inputs }

        //
        // Combine the genotyped VCFs from each family back together
        //

        VCF_GATHER_BCFTOOLS(
            ch_gather_inputs,
            ch_gather_inputs.map { meta, vcf, tbi ->
                [ meta, [], meta.region_count ]
            },
            "family",
            true
        )
        ch_versions = ch_versions.mix(VCF_GATHER_BCFTOOLS.out.versions)

        VCF_GATHER_BCFTOOLS.out.vcf
            .map { meta, vcf ->
                // Remove the bed counter from the meta field
                new_meta = meta - meta.subMap("region_count")
                [ new_meta, vcf ]
            }
            .set { ch_genotyped_vcfs }
    }

    emit:
    genotyped_vcfs = ch_genotyped_vcfs  // [ val(meta), path(vcf) ]
    versions       = ch_versions        // [ path(versions) ]
}
