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
    // Split BED file into correct regions
    //

    BEDTOOLS_SPLIT(
        beds.map { it + [params.scatter_count]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    BEDTOOLS_SPLIT.out.beds
        .set { split_beds }

    //
    // Generate DRAGSTR models
    //

    if (params.use_dragstr_model) {

        CALIBRATEDRAGSTRMODEL(
            crams,
            fasta,
            fasta_fai,
            dict,
            strtablefile
        )

        ch_versions = ch_versions.mix(CALIBRATEDRAGSTRMODEL.out.versions)

        crams
            .join(split_beds, failOnDuplicate: true, failOnMismatch: true)
            .join(CALIBRATEDRAGSTRMODEL.out.dragstr_model, failOnDuplicate: true, failOnMismatch: true)
            .set { cram_models }
    }
    else {
        crams
            .join(split_beds, failOnDuplicate: true, failOnMismatch: true)
            .set { cram_models }
    }

    //
    // Remap CRAM channel to fit the haplotypecaller input format
    //

    cram_models
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
        .map { meta, vcf, tbi ->
            new_meta = meta - meta.subMap("region")
            [ new_meta, vcf, tbi ]
        }
        .dump(tag:'haplotypecaller_vcfs', pretty:true)
        .set { haplotypecaller_vcfs }
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

    //
    // Merge the GVCFs
    //

    VCF_GATHER_BCFTOOLS(
        haplotypecaller_vcfs,
        haplotypecaller_vcfs.map { meta, vcf, tbi ->
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
        .tap { stats_input }
        .dump(tag:'gvcfs', pretty: true)
        .set { gvcfs }
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

    emit:
    gvcfs
    versions = ch_versions
    reports  = ch_reports
}
