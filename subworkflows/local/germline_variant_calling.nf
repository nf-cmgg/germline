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
    // Split the BED files into multiple subsets
    //

    beds
        .tap { dragstrmodel_beds }
        .set { beds_to_scatter }

    BED_SCATTER_GROOVY(
        beds_to_scatter.map {meta, bed -> [meta, bed]},
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
        crams
            .join(dragstrmodel_beds, failOnDuplicate: true, failOnMismatch: true)
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

        crams
            .combine(haplotypecaller_beds, by:0)
            .combine(CALIBRATEDRAGSTRMODEL.out.dragstr_model, by:0)
            .set { cram_models }
    }
    else {
        crams
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
        .join(HAPLOTYPECALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
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
        .join(VCF_GATHER_BCFTOOLS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
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
