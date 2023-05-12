//
// GERMLINE VARIANT CALLING
//

include { BEDTOOLS_SPLIT                                        } from '../../modules/nf-core/bedtools/split/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER              } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_CALIBRATEDRAGSTRMODEL as CALIBRATEDRAGSTRMODEL  } from '../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { DEEPVARIANT                                           } from '../../modules/nf-core/deepvariant/main'
include { TABIX_TABIX                                           } from '../../modules/nf-core/tabix/tabix/main'
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

    def callers  = params.callers.tokenize(",")
    ch_versions  = Channel.empty()
    ch_reports   = Channel.empty()
    ch_gvcfs     = Channel.empty()

    if("haplotypecaller" in callers) {
        //
        // Split BED file into multiple BEDs specified by --scatter_count
        //

        BEDTOOLS_SPLIT(
            ch_beds.map { it + [params.scatter_count] }
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

        BEDTOOLS_SPLIT.out.beds
            .set { ch_split_beds }

        //
        // Generate DRAGSTR models (if --dragstr true is specified)
        //

        if (params.dragstr) {

            CALIBRATEDRAGSTRMODEL(
                ch_crams,
                ch_fasta,
                ch_fai,
                ch_dict,
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
                new_meta = meta + [region_count: bed_is_list ? beds.size() : 1, caller:"haplotypecaller"]
                [ new_meta, cram, crai, bed_is_list ? beds : [beds], dragstr_model ]
            }
            .transpose(by:3) // Create one channel entry for each BED file per sample
            .map { meta, cram, crai, bed, dragstr_model ->
                // Set the base name of the BED file as the ID (this will look like sample_id.xxxx, where xxxx are numbers)
                new_meta = meta + [id:bed.baseName]
                [ new_meta, cram, crai, bed, dragstr_model ]
            }
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
        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions.first())

        HAPLOTYPECALLER.out.vcf
            .join(HAPLOTYPECALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .dump(tag:'haplotypecaller_vcfs', pretty:true)
            .set { ch_haplotypecaller_vcfs }

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
        ch_versions = ch_versions.mix(VCF_GATHER_BCFTOOLS.out.versions)

        VCF_GATHER_BCFTOOLS.out.vcf
            .join(VCF_GATHER_BCFTOOLS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi ->
                // Remove the bed counter from the meta field
                new_meta = meta - meta.subMap("region_count")
                [ new_meta, vcf, tbi ]
            }
            .dump(tag:'gvcfs', pretty: true)
            .set { ch_haplotypecaller_ready }

        ch_gvcfs = ch_gvcfs.mix(ch_haplotypecaller_ready)
    } 
    
    if("deepvariant" in callers) {
        ch_crams
            .join(ch_beds, failOnDuplicate:true, failOnMismatch:true)
            .map { meta, cram, crai, bed ->
                new_meta = meta + [caller:"deepvariant"]
                [ new_meta, cram, crai, bed ]
            }
            .set { ch_deepvariant_input }

        DEEPVARIANT(
            ch_deepvariant_input,
            ch_fasta,
            ch_fai,
            []
        )
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

        TABIX_TABIX(
            DEEPVARIANT.out.gvcf
        )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        DEEPVARIANT.out.gvcf
            .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .set { ch_deepvariant_ready }

        ch_gvcfs = ch_gvcfs.mix(ch_deepvariant_ready)
    }

    //
    // Run statistics on each GVCF produced
    //

    BCFTOOLS_STATS_INDIVIDUALS(
        ch_gvcfs,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS_INDIVIDUALS.out.versions.first())
    ch_reports  = ch_reports.mix(BCFTOOLS_STATS_INDIVIDUALS.out.stats.collect{it[1]})

    emit:
    gvcfs    = ch_gvcfs     // [ val(meta), path(gvcf), path(tbi) ]
    versions = ch_versions  // [ path(versions) ]
    reports  = ch_reports   // [ path(reports) ]
}
