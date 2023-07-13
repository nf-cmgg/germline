//
// GERMLINE VARIANT CALLING
//

include { BEDTOOLS_SPLIT                                        } from '../../modules/nf-core/bedtools/split/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER              } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { SAMTOOLS_CONVERT                                      } from '../../modules/nf-core/samtools/convert/main'
include { VARDICTJAVA                                           } from '../../modules/nf-core/vardictjava/main'
include { TABIX_TABIX                                           } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VCFS                             } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX                                      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_CONCAT                                       } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_REHEADER                                     } from '../../modules/nf-core/bcftools/reheader/main'
include { GATK4_CALIBRATEDRAGSTRMODEL as CALIBRATEDRAGSTRMODEL  } from '../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INDIVIDUALS          } from '../../modules/nf-core/bcftools/stats/main'

include { VCF_GATHER_BCFTOOLS                                   } from '../../subworkflows/nf-core/vcf_gather_bcftools/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their indexes
        ch_beds              // channel: [mandatory] [ val(meta), path(bed) ] => bed files created with the sample preparation subworkflow
        ch_fasta             // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_dict              // channel: [mandatory] [ path(dict) ] => sequence dictionary
        ch_strtablefile      // channel: [optional]  [ path(strtablefile) ] => STR table file
        ch_dbsnp             // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi         // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF

    main:

    ch_versions  = Channel.empty()
    ch_reports   = Channel.empty()

    val_callers = params.callers.tokenize(",")

    //
    // Split BED file into multiple BEDs specified by --scatter_count
    //

    BEDTOOLS_SPLIT(
        ch_beds.map { it + [params.scatter_count]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    BEDTOOLS_SPLIT.out.beds
        .set { ch_split_beds }

    ch_called_variants = Channel.empty()

    if("haplotypecaller" in val_callers) {

        //
        // Generate DRAGSTR models (if --dragstr is specified)
        //

        if (params.dragstr) {

            CALIBRATEDRAGSTRMODEL(
                ch_crams,
                ch_fasta.map { it[1] },
                ch_fai.map { it[1] },
                ch_dict.map { it[1] },
                ch_strtablefile
            )
            ch_versions = ch_versions.mix(CALIBRATEDRAGSTRMODEL.out.versions.first())

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
                // Determine the amount of BED files per sample
                bed_is_list = beds instanceof ArrayList
                new_meta = meta + [region_count: bed_is_list ? beds.size() : 1]
                [ new_meta, cram, crai, bed_is_list ? beds : [beds], dragstr_model ]
            }
            .transpose(by:3) // Create one channel entry for each BED file per sample
            .map { meta, cram, crai, bed, dragstr_model ->
                // Set the base name of the BED file as the ID (this will look like sample_id.xxxx, where xxxx are numbers)
                new_meta = meta + [id:bed.baseName, caller:"haplotypecaller"]
                [ new_meta, cram, crai, bed, dragstr_model ]
            }
            .dump(tag:'haplotypecaller_input', pretty:true)
            .set { ch_haplotypecaller_input }

        //
        // Call the variants using HaplotypeCaller
        //

        HAPLOTYPECALLER(
            ch_haplotypecaller_input,
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] },
            ch_dict.map { it[1] },
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions.first())

        ch_called_variants.mix(
            HAPLOTYPECALLER.out.vcf
                .join(HAPLOTYPECALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
                .dump(tag:'haplotypecaller_vcfs', pretty:true)
            )
            .set { ch_called_variants }
    }

    if("vardict" in val_callers) {
        //
        // Convert CRAMs to BAMs
        //

        ch_crams
            .branch { meta, cram, crai ->
                cram: cram.extension == "cram"
                bam: cram.extension == "bam"
            }
            .set { ch_cram_bam_branch }

        SAMTOOLS_CONVERT(
            ch_cram_bam_branch.cram,
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] }
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

        ch_cram_bam_branch.bam
            .mix(SAMTOOLS_CONVERT.out.alignment_index)
            .join(ch_split_beds, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, bam, bai, beds ->
                // Determine the amount of BED files per sample
                bed_is_list = beds instanceof ArrayList
                new_meta = meta + [region_count: bed_is_list ? beds.size() : 1]
                [ new_meta, bam, bai, bed_is_list ? beds : [beds] ]
            }
            .transpose(by:3) // Create one channel entry for each BED file per sample
            .map { meta, bam, bai, bed ->
                // Set the base name of the BED file as the ID (this will look like sample_id.xxxx, where xxxx are numbers)
                new_meta = meta + [id:bed.baseName, caller:"vardict"]
                [ new_meta, bam, bai, bed ]
            }
            .dump(tag:'vardict_input', pretty:true)
            .set { ch_vardict_input }

        //
        // Call the variants using Vardict
        //

        VARDICTJAVA(
            ch_vardict_input,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

        TABIX_BGZIPTABIX(
            VARDICTJAVA.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

        ch_called_variants.mix(TABIX_BGZIPTABIX.out.gz_tbi)
            .set { ch_called_variants }
    }

    //
    // Merge the GVCFs/VCFs produced by each caller
    //

    ch_called_variants
        .map { meta, vcf, tbi ->
            new_meta = meta + [id:meta.sample]
            [ groupKey(new_meta, meta.region_count), vcf, tbi ]
        }
        .groupTuple()
        .set { ch_concat_input }

    BCFTOOLS_CONCAT(
        ch_concat_input,
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_CONCAT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi ->
            // Remove the bed counter from the meta field
            new_meta = meta - meta.subMap("region_count")
            [ new_meta, vcf, tbi ]
        }
        .tap { ch_all_vcfs }
        .branch { meta, vcf, tbi ->
            gvcf: vcf.toString().endsWith("g.vcf.gz")
            vcf:  vcf.toString().endsWith("vcf.gz")
                return [ meta, vcf, file("${projectDir}/assets/vardict.header.vcf.gz", checkIfExists:true), tbi ]
        }
        .set { ch_vcfs }

    BCFTOOLS_REHEADER(
        ch_vcfs.vcf,
        ch_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions.first())

    TABIX_VCFS(
        BCFTOOLS_REHEADER.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_VCFS.out.versions.first())

    BCFTOOLS_REHEADER.out.vcf
        .join(TABIX_VCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .set { ch_vcfs_ready }

    //
    // Run statistics on each GVCF/VCF produced
    //

    BCFTOOLS_STATS_INDIVIDUALS(
        ch_vcfs.gvcf.mix(ch_vcfs_ready),
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS_INDIVIDUALS.out.versions.first())
    ch_reports  = ch_reports.mix(BCFTOOLS_STATS_INDIVIDUALS.out.stats.collect{it[1]})

    emit:
    gvcfs    = ch_vcfs.gvcf  // [ val(meta), path(gvcf), path(tbi) ]
    vcfs     = ch_vcfs_ready // [ val(meta), path(vcf), path(tbi) ]
    versions = ch_versions   // [ path(versions) ]
    reports  = ch_reports    // [ path(reports) ]
}
