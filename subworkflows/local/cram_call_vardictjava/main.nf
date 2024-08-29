include { SAMTOOLS_CONVERT                  } from '../../../modules/nf-core/samtools/convert/main'
include { VARDICTJAVA                       } from '../../../modules/nf-core/vardictjava/main'
include { TABIX_BGZIPTABIX as TABIX_SPLIT   } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_VCFANNO } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_REHEADER                 } from '../../../modules/nf-core/bcftools/reheader/main'
include { VCFANNO                           } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_TABIX                       } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS                    } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_CONCAT_BCFTOOLS               } from '../vcf_concat_bcftools/main'
include { VCF_FILTER_BCFTOOLS               } from '../vcf_filter_bcftools/main'

workflow CRAM_CALL_VARDICTJAVA {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their indexes
        ch_input             // channel: [mandatory] [ val(meta), path(cram), path(crai), path(bed) ] => sample CRAM files and their indexes
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dbsnp             // channel: [optional]  [ path(vcf) ] => the dbnsp vcf file
        ch_dbsnp_tbi         // channel: [optional]  [ path(tbi) ] => the dbsnp vcf index file
        filter               // boolean: filter the VCFs

    main:
        ch_versions = Channel.empty()

        ch_crams
            .map { meta, cram, crai ->
                def new_meta = meta + [caller:"vardict"]
                [ new_meta, cram, crai ]
            }
            .set { ch_crams }

        ch_crams
            .branch { meta, cram, crai ->
                bam: cram.extension == "bam"
                cram: cram.extension == "cram"
            }
            .set { ch_cram_bam }

        SAMTOOLS_CONVERT(
            ch_cram_bam.cram,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

        ch_input
            .map { meta, cram, crai, bed ->
                def new_meta = meta - meta.subMap("split_count") + [caller:"vardict", id:meta.sample]
                [ new_meta, cram, crai, bed, meta.split_count ]
            }
            .set { ch_vardict_crams }

        ch_cram_bam.bam
            .mix(SAMTOOLS_CONVERT.out.alignment_index)
            .combine(ch_vardict_crams, by:0)
            .map { meta, bam, bai, cram, crai, bed, split_count ->
                def new_meta = meta + [id:bed.baseName, split_count:split_count]
                [ new_meta, bam, bai, bed ]
            }
            .set { ch_vardict_input }

        VARDICTJAVA(
            ch_vardict_input,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

        TABIX_SPLIT(
            VARDICTJAVA.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_SPLIT.out.versions.first())

        VCF_CONCAT_BCFTOOLS(
            TABIX_SPLIT.out.gz_tbi,
            false
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

        ch_dbsnp_annotated = Channel.empty()
        if(ch_dbsnp) {
            ch_dbsnp
                .map { meta, dbsnp -> [ get_vcfanno_config(dbsnp) ] }
                .collect()
                .set { ch_vcfanno_toml }

            ch_dbsnp
                .combine(ch_dbsnp_tbi)
                .collect()
                .set { ch_vcfanno_resources }

            VCFANNO(
                VCF_CONCAT_BCFTOOLS.out.vcfs.map { meta, vcf -> [ meta, vcf, [], [] ] },
                ch_vcfanno_toml,
                [],
                ch_vcfanno_resources
            )
            ch_versions = ch_versions.mix(VCFANNO.out.versions.first())

            TABIX_VCFANNO(
                VCFANNO.out.vcf
            )
            ch_versions = ch_versions.mix(TABIX_VCFANNO.out.versions.first())

            TABIX_VCFANNO.out.gz_tbi.set { ch_dbsnp_annotated }
        } else {
            VCF_CONCAT_BCFTOOLS.out.vcfs.set { ch_dbsnp_annotated }
        }

        if(filter) {
            VCF_FILTER_BCFTOOLS(
                ch_dbsnp_annotated,
                false
            )
            ch_versions = ch_versions.mix(VCF_FILTER_BCFTOOLS.out.versions)
            ch_filter_output = VCF_FILTER_BCFTOOLS.out.vcfs
        } else {
            ch_filter_output = ch_dbsnp_annotated
        }

        TABIX_TABIX(
            ch_filter_output
        )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        ch_filter_output
            .join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi ->
                def new_meta = meta + [samples: meta.sample]
                [ new_meta, vcf, tbi ]
            }
            .set { ch_vcfs }

    emit:
    vcfs = ch_vcfs          // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions  // channel: [ path(versions.yml) ]

}

def get_vcfanno_config(vcf) {
    def old_toml = file("${projectDir}/assets/dbsnp.toml", checkIfExists: true)
    old_toml.copyTo("${workDir}/vcfanno/dbsnp.toml")
    def new_toml = file("${workDir}/vcfanno/dbsnp.toml")
    new_toml.text = old_toml.text.replace("DBSNP_FILE", vcf.getName())
    return new_toml
}
