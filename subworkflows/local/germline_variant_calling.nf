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
    // Generate DRAGSTR models
    //

    if (params.use_dragstr_model) {
        crams
            .join(beds, failOnDuplicate: true, failOnMismatch: true)
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
            .join(beds, failOnDuplicate: true, failOnMismatch: true)
            .join(CALIBRATEDRAGSTRMODEL.out.dragstr_model, failOnDuplicate: true, failOnMismatch: true)
            .set { cram_models }
    }
    else {
        crams
            .join(beds, failOnDuplicate: true, failOnMismatch: true)
            .set { cram_models }
    }

    //
    // Remap CRAM channel to fit the haplotypecaller input format
    //

    // TODO improve ROI handling => look at bedtools intersect being able to remove reads from bam/cram
    cram_models
        .dump(tag:'cram_models', pretty:true)
        .splitText(elem:3)
        .map(
            { meta, cram, crai, regions, dragstr_model=[] ->
                new_meta = meta + [region_count:regions.size()]
                [ new_meta, cram, crai, regions, dragstr_model ]
            }
        )
        .transpose()
        .map(
            { meta, cram, crai, region, dragstr_model ->
                region_split = region[0..-2].split("\t")
                region_string = "${region_split[0]}:${region_split[1]}-${region_split[2] as int + 1}"
                new_meta = meta + [id:"${meta.id}-${region_string.replace(":","_").replace("-":"_")}", region:region_string]
                [ new_meta, cram, crai, [], dragstr_model ]
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
        haplotypecaller_vcfs.map { meta, vcf ->
            [ meta, meta.region_count]
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
