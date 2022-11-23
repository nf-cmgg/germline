//
// GENOTYPE
//

include { RTGTOOLS_PEDFILTER as PEDFILTER            } from '../../modules/local/rtgtools/pedfilter/main'
include { MERGE_VCF_HEADERS                          } from '../../modules/local/merge_vcf_headers'
include { MERGE_BEDS                                 } from '../../modules/local/merge_beds'

include { GATK4_GENOMICSDBIMPORT as GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS as GENOTYPE_GVCFS      } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { TABIX_TABIX as TABIX_COMBINED_GVCFS        } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_POSTPROCESSED_VCFS    } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP as BGZIP_GENOTYPED_VCFS        } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_PED_VCFS              } from '../../modules/nf-core/tabix/bgzip/main'
include { BCFTOOLS_FILTER as FILTER_SNPS             } from '../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_INDELS           } from '../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_CONVERT                           } from '../../modules/nf-core/bcftools/convert/main'
include { BCFTOOLS_VIEW                              } from '../../modules/nf-core/bcftools/view/main'

workflow JOINT_GENOTYPING {
    take:
        gvcfs               // channel: [mandatory] [ meta, gvcf, tbi ] => The fresh GVCFs called with HaplotypeCaller
        beds                // channel: [mandatory] [ meta, bed ] => The BED files of the individuals
        peds                // channel: [mandatory] [ meta, peds ] => The pedigree files for the samples
        fasta               // channel: [mandatory] [ fasta ] => fasta reference
        fasta_fai           // channel: [mandatory] [ fasta_fai ] => fasta reference index
        dict                // channel: [mandatory] [ dict ] => sequence dictionary
        output_mode         // value:   [mandatory] whether or not to make the output seqplorer- or seqr-compatible

    main:

    ch_versions         = Channel.empty()

    beds
        .map(
            { meta, bed ->
                new_meta = [:]
                new_meta.family = meta.family
                new_meta.id = meta.family ?: meta.sample
                new_meta.family_count = meta.family_count
                [ groupKey(new_meta, meta.family_count.toInteger()), bed ]
            }
        )
        .groupTuple()
        .dump(tag:'bedmerge_input', pretty: true)
        .set { bedmerge_input }

    MERGE_BEDS(
        bedmerge_input
    )

    ch_versions = ch_versions.mix(MERGE_BEDS.out.versions)

    gvcfs
        .map(
            { meta, gvcf, tbi ->
                new_meta = [:]
                new_meta.family = meta.family
                new_meta.id = meta.family ?: meta.sample
                new_meta.family_count = meta.family_count
                [ groupKey(new_meta, meta.family_count.toInteger()), gvcf, tbi ]
            }
        )
        .groupTuple()
        .join(MERGE_BEDS.out.bed)
        .map(
            { meta, gvcfs, tbis, bed ->
                [ meta, gvcfs, tbis, bed, [], [] ]
            }
        )
        .set { genomicsdbimport_input }

    genomicsdbimport_input.dump(tag:'genomicsdbimport_input', pretty:true)

    //
    // Merge/Combine all the GVCFs from each family
    //

    GENOMICSDBIMPORT(
        genomicsdbimport_input,
        false,
        false,
        false
    )

    ch_versions = ch_versions.mix(GENOMICSDBIMPORT.out.versions)

    // TODO: Add BED splitting
    GENOMICSDBIMPORT.out.genomicsdb
        .map(
            { meta, db ->
                [ meta, db, [], [], [] ]
            }
        )
        .dump(tag:'genotypegvcfs_input', pretty:true)
        .set { genotypegvcfs_input }

    //
    // Genotype the combined GVCFs
    //

    GENOTYPE_GVCFS(
        genotypegvcfs_input,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    GENOTYPE_GVCFS.out.vcf.set { converted_vcfs }
    ch_versions = ch_versions.mix(GENOTYPE_GVCFS.out.versions)

    //
    // Add pedigree information
    //

    converted_vcfs
        .join(peds)
        .branch(
            { meta, vcf, ped ->
                has_ped: ped
                    return [ meta, vcf, ped ]
                no_ped: !ped
                    return [ meta, vcf ]
            }
        )
        .set { ped_vcfs }

    ped_vcfs.has_ped
        .tap { pedfilter_input }
        .dump(tag:'ped_vcfs_has_ped', pretty:true)
        .set { merge_header_input }
    ped_vcfs.no_ped.dump(tag:'ped_vcfs_no_ped', pretty:true)

    PEDFILTER(
        pedfilter_input.map({ meta, vcf, ped -> [ meta, ped ]})
    )

    ch_versions = ch_versions.mix(PEDFILTER.out.versions)

    BGZIP_GENOTYPED_VCFS(
        merge_header_input.map({ meta, vcf, ped -> [ meta, vcf ]})
    )

    ch_versions = ch_versions.mix(BGZIP_GENOTYPED_VCFS.out.versions)

    MERGE_VCF_HEADERS(
        BGZIP_GENOTYPED_VCFS.out.output
            .join(PEDFILTER.out.vcf)
    )

    ch_versions = ch_versions.mix(MERGE_VCF_HEADERS.out.versions)

    BGZIP_PED_VCFS(
        MERGE_VCF_HEADERS.out.vcf
    )

    ch_versions = ch_versions.mix(BGZIP_PED_VCFS.out.versions)

    ped_vcfs.no_ped
        .mix(BGZIP_PED_VCFS.out.output)
        .dump(tag:'filter_input', pretty:true)
        .set { filter_input }

    //
    // Filter the variants
    //

    // TODO: Move filter steps

    if (output_mode == "seqplorer") {
        FILTER_SNPS(
            filter_input
        )

        FILTER_INDELS(
            FILTER_SNPS.out.vcf
        )

        ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)
        ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

        FILTER_INDELS.out.vcf.set { genotyped_vcfs }

    }
    else {
        filter_input.set { genotyped_vcfs }
    }

    emit:
    genotyped_vcfs     // channel: [meta, vcf] => The output channel containing the post processed VCF
    versions = ch_versions
}
