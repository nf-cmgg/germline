//
// GERMLINE VARIANT CALLING
//

include { MERGE_BEDS                                            } from '../../modules/local/merge_beds'
include { SAMTOOLS_MERGE                                        } from '../../modules/local/samtools_merge'

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER              } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_CALIBRATEDRAGSTRMODEL as CALIBRATEDRAGSTRMODEL  } from '../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_REBLOCKGVCF as REBLOCKGVCF                      } from '../../modules/nf-core/gatk4/reblockgvcf/main'
include { BCFTOOLS_CONCAT                                       } from '../../modules/nf-core/bcftools/concat/main'
include { SAMTOOLS_INDEX                                        } from '../../modules/nf-core/samtools/index/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INDIVIDUALS          } from '../../modules/nf-core/bcftools/stats/main'
include { TABIX_TABIX as TABIX_GVCFS                            } from '../../modules/nf-core/tabix/tabix/main'

include { BED_SCATTER_BEDTOOLS                                  } from '../../subworkflows/nf-core/bed_scatter_bedtools/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        crams             // channel: [mandatory] [ meta, cram, crai ] => sample CRAM files and their indexes
        beds              // channel: [mandatory] [ meta, bed ] => bed files
        fasta             // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai         // channel: [mandatory] [ fasta_fai ] => fasta reference index
        dict              // channel: [mandatory] [ dict ] => sequence dictionary
        strtablefile      // channel: [mandatory] [ strtablefile ] => STR table file
        scatter_count     // value:   [mandatory] how many times the BED files need to be split before the variant calling
        use_dragstr_model // boolean: [mandatory] whether or not to use the dragstr models for variant calling
        always_use_cram   // boolean: [mandatory] whether or not to retain the bam after merging or convert back to cram
        dbsnp             // channel: [optional] The VCF containing the dbsnp variants
        dbsnp_tbi         // channel: [optional] The index of the dbsnp VCF

    main:

    gvcfs        = Channel.empty()
    ch_versions  = Channel.empty()
    ch_reports   = Channel.empty()

    //
    // Merge the CRAM files if there are multiple per sample
    //

    crams
        .groupTuple()
        .branch(
            { meta, cram, crai ->
                multiple: cram.size() > 1
                    return [meta, cram]
                single:   cram.size() == 1
                    return [meta, cram, crai]
            }
        )
        .set { cram_branch }

    cram_branch.multiple.dump(tag:'cram_branch_multiple', pretty:true)
    cram_branch.single.dump(tag:'cram_branch_single', pretty:true)

    SAMTOOLS_MERGE(
        cram_branch.multiple,
        fasta,
        fasta_fai,
        always_use_cram
    )

    SAMTOOLS_MERGE.out.cram
        .mix(SAMTOOLS_MERGE.out.bam)
        .mix(cram_branch.single
            .map(
                {meta, cram, crai ->
                    [ meta, cram[0], crai[0]]
                }
            )
        )
        .branch(
            { meta, cram, crai=[] ->
                not_indexed: crai == []
                    return [ meta, cram ]
                indexed: crai != []
                    return [ meta, cram, crai ]
            }
        )
        .set { merged_crams }

    merged_crams.not_indexed.dump(tag:'merged_crams_not_indexed', pretty:true)
    merged_crams.indexed.dump(tag:'merged_crams_indexed', pretty:true)

    SAMTOOLS_INDEX(
        merged_crams.not_indexed
    )

    merged_crams.not_indexed
        .join(SAMTOOLS_INDEX.out.crai)
        .mix(merged_crams.not_indexed
            .join(SAMTOOLS_INDEX.out.bai)
        )
        .mix(merged_crams.indexed)
        .dump(tag:'ready_crams', pretty:true)
        .set { ready_crams }

    //
    // Merge the BED files if there are multiple per sample
    //

    beds
        .groupTuple()
        .branch(
            { meta, bed ->
                multiple: bed.size() > 1
                    return [meta, bed]
                single:   bed.size() == 1
                    return [meta, bed]
            }
        )
        .set { bed_branch }

    MERGE_BEDS(
        bed_branch.multiple
    )

    MERGE_BEDS.out.bed
        .mix(bed_branch.single)
        .map(
            { meta, bed ->
                [ meta, bed, scatter_count ]
            }
        )
        .dump(tag:'merged_beds', pretty:true)
        .set { merged_beds }

    //
    // Split the BED files into multiple subsets
    //

    if (scatter_count > 1) {
        BED_SCATTER_BEDTOOLS(
            merged_beds
        )

        ch_versions = ch_versions.mix(BED_SCATTER_BEDTOOLS.out.versions)

        BED_SCATTER_BEDTOOLS.out.scattered_beds
            .set { split_beds }
    }
    else {
        merged_beds
        .map(
            { meta, bed, bed_count ->
                [ meta, bed, 1 ]
            }
        )
        .set { split_beds }
    }

    split_beds.dump(tag:'split_beds', pretty:true)

    //
    // Generate DRAGSTR models
    //

    if (use_dragstr_model) {
        ready_crams
            .join(merged_beds)
            .map(
                { meta, cram, crai, bed, bed_count ->
                    [meta, cram, crai, bed]
                }
            )
            .set { calibratedragstrmodel_input }

        CALIBRATEDRAGSTRMODEL(
            calibratedragstrmodel_input,
            fasta,
            fasta_fai,
            dict,
            strtablefile
        )

        ch_versions = ch_versions.mix(CALIBRATEDRAGSTRMODEL.out.versions)

        ready_crams
            .combine(split_beds, by:0)
            .combine(CALIBRATEDRAGSTRMODEL.out.dragstr_model, by:0)
            .set { cram_models }
    }
    else {
        ready_crams
            .combine(split_beds, by:0)
            .set { cram_models }
    }

    cram_models.dump(tag:'cram_models', pretty:true)

    //
    // Remap CRAM channel to fit the haplotypecaller input format
    //

    cram_models
        .map(
            { meta, cram, crai, bed, bed_count, dragstr_model=[] ->
                new_meta = meta.clone()

                // If either no scatter is done, i.e. one interval (1), then don't rename samples
                new_meta.id = bed_count <= 1 ? meta.id : bed.baseName
                new_meta.bed_count = bed_count

                [ new_meta, cram, crai, bed, dragstr_model ]
            }
        )
        .dump(tag:'cram_intervals', pretty:true)
        .set { cram_intervals }

    //
    // Call the variants using HaplotypeCaller
    //

    HAPLOTYPECALLER(
        cram_intervals,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    HAPLOTYPECALLER.out.vcf
        .join(HAPLOTYPECALLER.out.tbi)
        .dump(tag:'haplotypecaller_vcfs', pretty:true)
        .set { haplotypecaller_vcfs }
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

    //
    // Merge the GVCFs if split BED files were used
    //

    if (scatter_count > 1) {
        haplotypecaller_vcfs
            .map(
                { meta, vcf, tbi ->
                    new_meta = meta.clone()
                    new_meta.id = new_meta.sample
                    new_meta.remove('bed_count')
                    [ groupKey(new_meta, meta.bed_count.toInteger()), vcf, tbi ]
                }
            )
            .groupTuple()
            .branch(
                { meta, vcfs, tbis ->
                    multiple: vcfs.size() > 1
                        return [ meta, vcfs, tbis ]
                    single: vcfs.size() == 1
                        return [ meta, vcfs[0] ]
                }
            )
            .set { concat_input }

        concat_input.multiple.dump(tag:'concat_input_multiple', pretty:true)
        concat_input.single.dump(tag:'concat_input_single', pretty:true)

        BCFTOOLS_CONCAT(
            concat_input.multiple
        )

        BCFTOOLS_CONCAT.out.vcf
            .mix(concat_input.single)
            .set { tabixgvcfs_input }
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
    }
    else {
        haplotypecaller_vcfs
            .map(
                { meta, vcf, tbi ->
                    [ meta, vcf ]
                }
            )
            .set { tabixgvcfs_input }
    }

    tabixgvcfs_input.dump(tag:'tabixgvcfs_input', pretty: true)

    //
    // Create indices for all the GVCF files
    //

    TABIX_GVCFS(
        tabixgvcfs_input
    )

    tabixgvcfs_input
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
    // Perform QC on the individual GVCF
    //

    BCFTOOLS_STATS_INDIVIDUALS(
        tabixgvcfs_input.join(TABIX_GVCFS.out.tbi),
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS_INDIVIDUALS.out.versions)
    ch_reports  = ch_reports.mix(BCFTOOLS_STATS_INDIVIDUALS.out.stats.collect{it[1]})

    //
    // Reblock the single sample GVCF files
    //

    REBLOCKGVCF(
        reblockgvcf_input,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    ch_versions = ch_versions.mix(REBLOCKGVCF.out.versions)

    emit:
    gvcfs    = REBLOCKGVCF.out.vcf
    versions = ch_versions
    reports  = ch_reports
}
