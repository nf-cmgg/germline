//
// Call the variants using HaplotypeCaller
//

include { GATK4_CALIBRATEDRAGSTRMODEL               } from '../../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_HAPLOTYPECALLER                     } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { BCFTOOLS_STATS                            } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_CONCAT_BCFTOOLS                       } from '../vcf_concat_bcftools/main'

workflow CRAM_CALL_HAPLOTYPECALLER {
    take:
        ch_input            // channel: [mandatory] [ val(meta), path(cram), path(crai), path(bed) ] => sample CRAM files and their indexes with the split bed files
        ch_fasta            // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai              // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dict             // channel: [mandatory] [ val(meta), path(dict) ] => sequence dictionary
        ch_strtablefile     // channel: [optional]  [ val(meta), path(strtablefile) ] => STR table file
        ch_dbsnp            // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi        // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF
        dragstr             // boolean: create a DragSTR model and run haplotypecaller with it

    main:

    //
    // Generate DRAGSTR models (if --dragstr is specified)
    //

    def ch_cram_models = channel.empty()
    if (dragstr) {

        ch_input
            .map { meta, cram, crai, bed ->
                def new_meta = meta + [id:meta.sample]
                [ groupKey(new_meta, meta.split_count), cram, crai, bed ]
            }
            .tap { ch_original }
            .groupTuple()
            .map { meta, cram, crai, beds ->
                [ meta, cram.unique()[0], crai.unique()[0], beds ]
            }
            .set { ch_dragstr_input }

        GATK4_CALIBRATEDRAGSTRMODEL(
            ch_dragstr_input.map { meta, cram, crai, _beds -> [ meta, cram, crai ] },
            ch_fasta.map { _meta, fasta -> fasta },
            ch_fai.map { _meta, fai -> fai },
            ch_dict.map { _meta, dict -> dict },
            ch_strtablefile.map { _meta, str -> str }
        )

        ch_cram_models = ch_original
            .combine(GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model, by: 0)
            .map { meta, cram, crai, bed, dragstr_model ->
                def new_meta = meta + [id:bed.baseName]
                [ new_meta, cram, crai, bed, dragstr_model ]
            }
    }
    else {
        ch_cram_models = ch_input
            .map { meta, cram, crai, bed -> [ meta, cram, crai, bed, [] ] }
    }

    GATK4_HAPLOTYPECALLER(
        ch_cram_models,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )

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

    BCFTOOLS_STATS(
        VCF_CONCAT_BCFTOOLS.out.vcfs,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )

    emit:
    gvcfs = VCF_CONCAT_BCFTOOLS.out.vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]
    reports = BCFTOOLS_STATS.out.stats      // channel: [ val(meta), path(stats) ]
}
