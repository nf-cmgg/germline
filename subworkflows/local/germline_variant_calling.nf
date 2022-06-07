//
// GERMLINE VARIANT CALLING
//

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER  } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { BCFTOOLS_CONCAT                           } from '../../modules/nf-core/modules/bcftools/concat/main'
include { BEDTOOLS_SPLIT                            } from '../../modules/nf-core/modules/bedtools/split/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        cram         // channel: [mandatory] cram
        beds         // channel: [mandatory] bed regions
        fasta        // channel: [mandatory] fasta reference
        fasta_fai    // channel: [mandatory] fasta reference index
        dict         // channel: [mandatory] sequence dictionary
        strtablefile // channel: [mandatory] STR table file

    main:

    vcfs             = Channel.empty()
    ch_versions      = Channel.empty()

    //
    // Split the BED files into multiple subsets
    //

    if(params.scatter_count > 1){

        BEDTOOLS_SPLIT(
            beds,
            params.scatter_count
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions)

        split_beds = BEDTOOLS_SPLIT.out.beds
        .transpose()
        .map({ meta, bed ->
            [ meta, bed, params.scatter_count ]
        })
    }
    else{
        split_beds = beds
        .map({ meta, bed ->
            [ meta, bed, params.scatter_count ]
        })
    }

    //
    // Remap channel with split_beds
    //

    cram_intervals = cram.combine(split_beds, by: 0)
        .map{ meta, cram, crai, split_beds, num_intervals ->
            new_meta = meta.clone()
            new_meta.sample = new_meta.id

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_meta.id = num_intervals <= 1 ? meta.id : split_beds.baseName
            new_meta.num_intervals = num_intervals

            //If no interval file provided (0) then add empty list
            split_beds_new = num_intervals == 0 ? [] : split_beds

            [new_meta, cram, crai, split_beds_new]
        }

    //
    // Call the variants using HaplotypeCaller
    //

    HAPLOTYPECALLER(
        cram_intervals,
        fasta,
        fasta_fai,
        dict,
        [],
        []
    )

    haplotypecaller_vcfs = HAPLOTYPECALLER.out.vcf.combine(HAPLOTYPECALLER.out.tbi, by:0)
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

    //
    // Merge the VCFs if split BED files were used
    //
    
    if (params.scatter_count > 1){

        concat_input = haplotypecaller_vcfs
                    .map({meta, vcf, tbi -> 
                        new_meta = meta.clone()
                        new_meta.id = new_meta.sample
                        [ new_meta, vcf, tbi ]
                    })
                    .groupTuple()

        BCFTOOLS_CONCAT(
            concat_input
        )

        vcfs = BCFTOOLS_CONCAT.out.vcf
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
    }
    else {
        vcfs = haplotypecaller_vcfs
                        .map({ meta, vcf, tbi ->
                            [ meta, vcf ]
                        })
    }

    emit:
    vcfs    
    versions = ch_versions
}