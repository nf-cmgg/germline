//
// GERMLINE VARIANT CALLING
//

include { BEDTOOLS_SPLIT                                        } from '../../modules/nf-core/bedtools/split/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER              } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_CALIBRATEDRAGSTRMODEL as CALIBRATEDRAGSTRMODEL  } from '../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INDIVIDUALS          } from '../../modules/nf-core/bcftools/stats/main'

include { VCF_GATHER_BCFTOOLS                                   } from '../../subworkflows/nf-core/vcf_gather_bcftools/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        ch_crams             // channel: [mandatory] [ meta, cram, crai ] => sample CRAM files and their indexes
        ch_beds              // channel: [mandatory] [ meta, bed ] => bed files
        ch_fasta             // channel: [mandatory] [ fasta ] => fasta reference
        ch_fai               // channel: [mandatory] [ fasta_fai ] => fasta reference index
        ch_dict              // channel: [mandatory] [ dict ] => sequence dictionary
        ch_strtablefile      // channel: [optional] [ strtablefile ] => STR table file
        ch_dbsnp             // channel: [optional] The VCF containing the dbsnp variants
        ch_dbsnp_tbi         // channel: [optional] The index of the dbsnp VCF

    main:

    ch_versions  = Channel.empty()
    ch_reports   = Channel.empty()

    //
    // Split BED file into correct regions
    //

    BEDTOOLS_SPLIT(
        ch_beds.map { it + [params.scatter_count]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    BEDTOOLS_SPLIT.out.beds
        .set { ch_split_beds }

    //
    // Generate DRAGSTR models
    //

    if (params.use_dragstr_model) {

        CALIBRATEDRAGSTRMODEL(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_strtablefile
        )

        ch_versions = ch_versions.mix(CALIBRATEDRAGSTRMODEL.out.versions)

        ch_crams
            .join(ch_split_beds, failOnDuplicate: true, failOnMismatch: true)
            .join(CALIBRATEDRAGSTRMODEL.out.dragstr_model, failOnDuplicate: true, failOnMismatch: true)
            .set { ch_cram_models }
    }
    else {
        ch_crams
            .join(ch_split_beds, failOnDuplicate: true, failOnMismatch: true)
            .set { ch_cram_models }
    }

    //
    // Remap CRAM channel to fit the haplotypecaller input format
    //

    ch_cram_models
        .dump(tag:'cram_models', pretty:true)
        .map { meta, cram, crai, beds, dragstr_model=[] ->
            bed_is_list = beds instanceof ArrayList
            new_meta = meta + [region_count: bed_is_list ? beds.size() : 1]
            [ new_meta, cram, crai, bed_is_list ? beds : [beds], dragstr_model ]
        }
        .transpose(by:3)
        .map(
            { meta, cram, crai, bed, dragstr_model ->
                new_meta = meta + [id:bed.baseName]
                [ new_meta, cram, crai, bed, dragstr_model ]
            }
        )
        .dump(tag:'cram_intervals', pretty:true)
        .set { ch_cram_intervals }

    //
    // Call the variants using HaplotypeCaller
    //

    HAPLOTYPECALLER(
        ch_cram_intervals,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )

    HAPLOTYPECALLER.out.vcf
        .join(HAPLOTYPECALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi ->
            new_meta = meta - meta.subMap("region")
            [ new_meta, vcf, tbi ]
        }
        .dump(tag:'haplotypecaller_vcfs', pretty:true)
        .set { ch_haplotypecaller_vcfs }
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

    //
    // Merge the GVCFs
    //

    VCF_GATHER_BCFTOOLS(
        ch_haplotypecaller_vcfs,
        ch_haplotypecaller_vcfs.map { meta, vcf, tbi ->
            [ meta, [], meta.region_count ]
        },
        "sample",
        false
    )

    VCF_GATHER_BCFTOOLS.out.vcf
        .join(VCF_GATHER_BCFTOOLS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi ->
            new_meta = meta - meta.subMap("region_count")
            [ new_meta, vcf, tbi ]
        }
        .dump(tag:'gvcfs', pretty: true)
        .set { ch_gvcfs }
    ch_versions = ch_versions.mix(VCF_GATHER_BCFTOOLS.out.versions)

    //
    // Create indices for all the GVCF files
    //

    BCFTOOLS_STATS_INDIVIDUALS(
        ch_gvcfs,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS_INDIVIDUALS.out.versions)
    ch_reports  = ch_reports.mix(BCFTOOLS_STATS_INDIVIDUALS.out.stats.collect{it[1]})

    emit:
    gvcfs    = ch_gvcfs
    versions = ch_versions
    reports  = ch_reports
}
