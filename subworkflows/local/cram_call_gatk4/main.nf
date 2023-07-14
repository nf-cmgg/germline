//
// Call the variants using GATK4 tooling
//

include { GATK4_CALIBRATEDRAGSTRMODEL   } from '../../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_HAPLOTYPECALLER         } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { BCFTOOLS_STATS                } from '../../../modules/nf-core/bcftools/stats/main'

include { INPUT_SPLIT_BEDTOOLS          } from '../input_split_bedtools/main'
include { VCF_CONCAT_BCFTOOLS           } from '../vcf_concat_bcftools/main'

workflow CRAM_CALL_GATK4 {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their indexes
        ch_beds              // channel: [mandatory] [ val(meta), path(bed) ] => bed files created with the sample preparation subworkflow
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dict              // channel: [mandatory] [ val(meta), path(dict) ] => sequence dictionary
        ch_strtablefile      // channel: [optional]  [ path(strtablefile) ] => STR table file
        ch_dbsnp             // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi         // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF

    main:

    ch_versions  = Channel.empty()

    //
    // Generate DRAGSTR models (if --dragstr is specified)
    //

    if (params.dragstr) {

        GATK4_CALIBRATEDRAGSTRMODEL(
            ch_crams,
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] },
            ch_dict.map { it[1] },
            ch_strtablefile
        )
        ch_versions = ch_versions.mix(GATK4_CALIBRATEDRAGSTRMODEL.out.versions.first())

        ch_crams
            .join(GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model, failOnDuplicate: true, failOnMismatch: true)
            .set { ch_cram_models }
    }
    else {
        ch_crams
            .map { meta, cram, crai -> [ meta, cram, crai, [] ] }
            .set { ch_cram_models }
    }

    INPUT_SPLIT_BEDTOOLS(
        ch_beds.map{ it + [params.scatter_count] },
        ch_cram_models
    )
    ch_versions = ch_versions.mix(INPUT_SPLIT_BEDTOOLS.out.versions)

    INPUT_SPLIT_BEDTOOLS.out.split
        .map { meta, cram, crai, dragstr, bed ->
            [ meta, cram, crai, bed, dragstr ]
        }
        .set { ch_split }

    GATK4_HAPLOTYPECALLER(
        ch_split,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] },
        ch_dict.map { it[1] },
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi ->
            def new_meta = meta + [caller: "haplotypecaller"]
            [ new_meta, vcf, tbi ]
        }
        .set { ch_called_variants }

    VCF_CONCAT_BCFTOOLS(
        ch_called_variants
    )
    ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

    BCFTOOLS_STATS(
        VCF_CONCAT_BCFTOOLS.out.vcfs,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
    gvcfs = VCF_CONCAT_BCFTOOLS.out.vcfs                // channel: [ val(meta), path(vcf), path(tbi) ]

    reports = BCFTOOLS_STATS.out.stats.collect{it[1]}   // channel: [ path(stats) ]
    versions = ch_versions                              // channel: [ versions.yml ]

}