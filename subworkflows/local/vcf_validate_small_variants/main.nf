include { RTGTOOLS_VCFEVAL      } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { RTGTOOLS_ROCPLOT      } from '../../../modules/nf-core/rtgtools/rocplot/main'

workflow VCF_VALIDATE_SMALL_VARIANTS {

    take:
    ch_vcf                          // [mandatory] channel: [ meta, vcf, tbi, truth_vcf, truth_tbi ]
    ch_beds                         // [mandatory] channel: [ meta, regions_bed, targets_bed ]
    ch_fasta                        // [happy only] channel: [ meta, fasta ]
    ch_fasta_fai                    // [happy only] channel: [ meta, fasta_fai ]
    ch_vcfeval_sdf                  // [vcfeval only] channel: [ meta, sdf ]

    main:

    ch_versions                             = Channel.empty()

    ch_input = ch_vcf.join(ch_beds, failOnDuplicate: true, failOnMismatch: true)

    RTGTOOLS_VCFEVAL(
        ch_input,
        ch_vcfeval_sdf
    )
    ch_versions = ch_versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

    ch_rocplot_input = RTGTOOLS_VCFEVAL.out.snp_roc
        .map { meta, tsv ->
            [ meta + [roc_type:'snp'], tsv ]
        }
        .mix(
            RTGTOOLS_VCFEVAL.out.non_snp_roc.map { meta, tsv ->
                [ meta + [roc_type:'non_snp'], tsv ]
            },
            RTGTOOLS_VCFEVAL.out.weighted_roc.map { meta, tsv ->
                [ meta + [roc_type:'weighted'], tsv ]
            }
        )

    vcfeval_true_positive_vcf               = RTGTOOLS_VCFEVAL.out.tp_vcf
    vcfeval_true_positive_vcf_tbi           = RTGTOOLS_VCFEVAL.out.tp_tbi
    vcfeval_false_negative_vcf              = RTGTOOLS_VCFEVAL.out.fn_vcf
    vcfeval_false_negative_vcf_tbi          = RTGTOOLS_VCFEVAL.out.fn_tbi
    vcfeval_false_positive_vcf              = RTGTOOLS_VCFEVAL.out.fp_vcf
    vcfeval_false_positive_vcf_tbi          = RTGTOOLS_VCFEVAL.out.fp_tbi
    vcfeval_true_positive_baseline_vcf      = RTGTOOLS_VCFEVAL.out.baseline_vcf
    vcfeval_true_positive_baseline_vcf_tbi  = RTGTOOLS_VCFEVAL.out.baseline_tbi
    vcfeval_summary                         = RTGTOOLS_VCFEVAL.out.summary
    vcfeval_phasing                         = RTGTOOLS_VCFEVAL.out.phasing
    vcfeval_snp_roc                         = RTGTOOLS_VCFEVAL.out.snp_roc
    vcfeval_non_snp_roc                     = RTGTOOLS_VCFEVAL.out.non_snp_roc
    vcfeval_weighted_roc                    = RTGTOOLS_VCFEVAL.out.weighted_roc

    RTGTOOLS_ROCPLOT(
        ch_rocplot_input
    )

    ch_versions = ch_versions.mix(RTGTOOLS_ROCPLOT.out.versions.first())

    rocplot_out_png = RTGTOOLS_ROCPLOT.out.png
        .branch { meta, png ->
            roc_type = meta.roc_type
            def new_meta = meta - meta.subMap("roc_type")

            snp:        roc_type == "snp"
            non_snp:    roc_type == "non_snp"
            weighted:   roc_type == "weighted"
        }

    rocplot_out_svg = RTGTOOLS_ROCPLOT.out.svg
        .branch { meta, svg ->
            roc_type = meta.roc_type
            def new_meta = meta - meta.subMap("roc_type")

            snp:        roc_type == "snp"
            non_snp:    roc_type == "non_snp"
            weighted:   roc_type == "weighted"
        }

    rtgtools_snp_png_rocplot        = rocplot_out_png.snp
    rtgtools_non_snp_png_rocplot    = rocplot_out_png.non_snp
    rtgtools_weighted_png_rocplot   = rocplot_out_png.weighted

    rtgtools_snp_svg_rocplot        = rocplot_out_svg.snp
    rtgtools_non_snp_svg_rocplot    = rocplot_out_svg.non_snp
    rtgtools_weighted_svg_rocplot   = rocplot_out_svg.weighted

    emit:
    vcfeval_true_positive_vcf               // channel: [ meta, vcf ]
    vcfeval_true_positive_vcf_tbi           // channel: [ meta, tbi ]
    vcfeval_false_negative_vcf              // channel: [ meta, vcf ]
    vcfeval_false_negative_vcf_tbi          // channel: [ meta, tbi ]
    vcfeval_false_positive_vcf              // channel: [ meta, vcf ]
    vcfeval_false_positive_vcf_tbi          // channel: [ meta, tbi ]
    vcfeval_true_positive_baseline_vcf      // channel: [ meta, vcf ]
    vcfeval_true_positive_baseline_vcf_tbi  // channel: [ meta, tbi ]
    vcfeval_summary                         // channel: [ meta, summary ]
    vcfeval_phasing                         // channel: [ meta, phasing ]
    vcfeval_snp_roc                         // channel: [ meta, tsv ]
    vcfeval_non_snp_roc                     // channel: [ meta, tsv ]
    vcfeval_weighted_roc                    // channel: [ meta, tsv ]

    rtgtools_snp_png_rocplot                // channel: [ meta, png ]
    rtgtools_non_snp_png_rocplot            // channel: [ meta, png ]
    rtgtools_weighted_png_rocplot           // channel: [ meta, png ]
    rtgtools_snp_svg_rocplot                // channel: [ meta, svg ]
    rtgtools_non_snp_svg_rocplot            // channel: [ meta, svg ]
    rtgtools_weighted_svg_rocplot           // channel: [ meta, svg ]

    versions = ch_versions                  // channel: [ versions.yml ]
}

