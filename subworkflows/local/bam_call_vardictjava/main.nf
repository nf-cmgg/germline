include { VARDICTJAVA                       } from '../../../modules/nf-core/vardictjava/main'
include { TABIX_BGZIP                       } from '../../../modules/nf-core/tabix/bgzip/main'
include { BCFTOOLS_REHEADER                 } from '../../../modules/nf-core/bcftools/reheader/main'
include { VCFANNO                           } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_TABIX                       } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS                    } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_CONCAT_BCFTOOLS               } from '../vcf_concat_bcftools/main'
include { VCF_FILTER_BCFTOOLS               } from '../vcf_filter_bcftools/main'

workflow BAM_CALL_VARDICTJAVA {
    take:
        ch_input             // channel: [mandatory] [ val(meta), path(bam), path(bai), path(bed) ] => sample CRAM files and their indexes
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dbsnp             // channel: [optional]  [ path(vcf) ] => the dbnsp vcf file
        ch_dbsnp_tbi         // channel: [optional]  [ path(tbi) ] => the dbsnp vcf index file
        filter               // boolean: filter the VCFs

    main:
    def ch_versions = Channel.empty()

    VARDICTJAVA(
        ch_input.map { meta, bam, bai, bed ->
            def new_meta = meta + [caller:'vardict']
            [ new_meta, bam, bai, bed ]
        },
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

    VCF_CONCAT_BCFTOOLS(
        VARDICTJAVA.out.vcf,
        false
    )
    ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

    def ch_annotated = Channel.empty()
    if(!(ch_dbsnp instanceof List)) {
        ch_dbsnp.map { _meta, dbsnp -> [ get_vcfanno_config(dbsnp) ] }
            .collect()
            .set { ch_vcfanno_toml } // Set needs to be used here due to some Nextflow bug

        ch_dbsnp.map { _meta, dbsnp -> dbsnp }
            .combine(ch_dbsnp_tbi.map { _meta, tbi -> tbi })
            .collect()
            .set { ch_vcfanno_resources } // Set needs to be used here due to some Nextflow bug

        VCFANNO(
            VCF_CONCAT_BCFTOOLS.out.vcfs.map { meta, vcf -> [ meta, vcf, [], [] ] },
            ch_vcfanno_toml,
            [],
            ch_vcfanno_resources
        )
        ch_versions = ch_versions.mix(VCFANNO.out.versions.first())

        TABIX_BGZIP(
            VCFANNO.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

        ch_annotated = TABIX_BGZIP.out.output
    } else {
        ch_annotated = VCF_CONCAT_BCFTOOLS.out.vcfs
    }

    def ch_filter_output = Channel.empty()
    if(filter) {
        VCF_FILTER_BCFTOOLS(
            ch_annotated,
            false
        )
        ch_versions = ch_versions.mix(VCF_FILTER_BCFTOOLS.out.versions)
        ch_filter_output = VCF_FILTER_BCFTOOLS.out.vcfs
    } else {
        ch_filter_output = ch_annotated
    }

    TABIX_TABIX(
        ch_filter_output
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    def ch_vcfs = ch_filter_output
        .join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi ->
            def new_meta = meta + [family_samples: meta.sample]
            [ new_meta, vcf, tbi ]
        }

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
