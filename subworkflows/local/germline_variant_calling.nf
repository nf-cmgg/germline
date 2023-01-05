//
// GERMLINE VARIANT CALLING
//

include { MERGE_BEDS                                            } from '../../modules/local/merge_beds'
include { SAMTOOLS_MERGE                                        } from '../../modules/local/samtools_merge'

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER              } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_CALIBRATEDRAGSTRMODEL as CALIBRATEDRAGSTRMODEL  } from '../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_REBLOCKGVCF as REBLOCKGVCF                      } from '../../modules/nf-core/gatk4/reblockgvcf/main'
include { SAMTOOLS_INDEX                                        } from '../../modules/nf-core/samtools/index/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INDIVIDUALS          } from '../../modules/nf-core/bcftools/stats/main'
include { TABIX_TABIX as TABIX_GVCFS                            } from '../../modules/nf-core/tabix/tabix/main'

include { BED_SCATTER_GROOVY                                    } from '../../subworkflows/local/bed_scatter_groovy'

include { VCF_GATHER_BCFTOOLS                                   } from '../../subworkflows/nf-core/vcf_gather_bcftools/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        crams             // channel: [mandatory] [ meta, cram, crai ] => sample CRAM files and their indexes
        beds              // channel: [mandatory] [ meta, bed ] => bed files
        fasta             // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai         // channel: [mandatory] [ fasta_fai ] => fasta reference index
        dict              // channel: [mandatory] [ dict ] => sequence dictionary
        strtablefile      // channel: [mandatory] [ strtablefile ] => STR table file
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
        fasta_fai
    )

    SAMTOOLS_MERGE.out.cram
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

    merged_crams.not_indexed
        .tap { crams_to_index }
        .dump(tag:'merged_crams_not_indexed', pretty:true)
        .set { crams_without_index }
    merged_crams.indexed.dump(tag:'merged_crams_indexed', pretty:true)

    SAMTOOLS_INDEX(
        crams_to_index
    )

    crams_without_index
        .join(SAMTOOLS_INDEX.out.crai)
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
        .tap { dragstrmodel_beds }
        .dump(tag:'merged_beds', pretty:true)
        .set { beds_to_scatter }

    //
    // Split the BED files into multiple subsets
    //

    BED_SCATTER_GROOVY(
        beds_to_scatter.map {meta, bed -> [meta, bed instanceof Path ? bed : bed[0]]},
        params.scatter_size
    )

    ch_versions = ch_versions.mix(BED_SCATTER_GROOVY.out.versions)

    BED_SCATTER_GROOVY.out.scattered
        .tap { gather_beds }
        .dump(tag:'split_beds', pretty:true)
        .set { haplotypecaller_beds }

    //
    // Generate DRAGSTR models
    //

    if (params.use_dragstr_model) {
        ready_crams
            .join(dragstrmodel_beds)
            .map(
                { meta, cram, crai, bed ->
                    [meta, cram, crai, bed]
                }
            )
            .dump(tag:'calibratedragstrmodel_input')
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
            .combine(haplotypecaller_beds, by:0)
            .combine(CALIBRATEDRAGSTRMODEL.out.dragstr_model, by:0)
            .set { cram_models }
    }
    else {
        ready_crams
            .combine(haplotypecaller_beds, by:0)
            .set { cram_models }
    }

    //
    // Remap CRAM channel to fit the haplotypecaller input format
    //

    cram_models
        .dump(tag:'cram_models', pretty:true)
        .map(
            { meta, cram, crai, bed, bed_count, dragstr_model=[] ->
                new_meta = meta.clone()
                new_meta.id = bed.baseName
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
    // Merge the GVCFs
    //

    VCF_GATHER_BCFTOOLS(
        haplotypecaller_vcfs,
        gather_beds.map { meta, bed, scatter_count ->
            new_meta = meta.findAll{true} + [id:bed.baseName]
            [ new_meta, bed, scatter_count ]
        },
        "sample",
        false
    )

    VCF_GATHER_BCFTOOLS.out.vcf
        .join(VCF_GATHER_BCFTOOLS.out.tbi)
        .tap { stats_input }
        .dump(tag:'reblockgvcf_input', pretty: true)
        .set { reblockgvcf_input }
    ch_versions = ch_versions.mix(VCF_GATHER_BCFTOOLS.out.versions)

    //
    // Create indices for all the GVCF files
    //

    BCFTOOLS_STATS_INDIVIDUALS(
        stats_input,
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
        reblockgvcf_input.map{ it + [[]] },
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
