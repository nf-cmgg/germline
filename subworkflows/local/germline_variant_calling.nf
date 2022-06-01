//
// GERMLINE VARIANT CALLING
//

include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'

workflow GERMLINE_VARIANT_CALLING {
    take:
        cram                         // channel: [mandatory] cram
        intervals                    // channel: [mandatory] intervals/target regions

    main:

    vcfs             = Channel.empty()
    ch_versions      = Channel.empty()

    // Remap channel with intervals
    cram_intervals = cram.combine(intervals, by: 0)
        .map{ meta, cram, crai, intervals, num_intervals ->
            new_meta = meta.clone()
            new_meta.sample = new_meta.id

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_meta.id = num_intervals <= 1 ? meta.id : intervals.baseName
            new_meta.num_intervals = num_intervals

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [new_meta, cram, crai, intervals_new]
        }

    HAPLOTYPECALLER(
        cram_intervals,
        params.fasta,
        params.fasta_fai,
        params.dict,
        [],
        []
    )

    haplotypecaller_vcfs = HAPLOTYPECALLER.out.vcf.combine(HAPLOTYPECALLER.out.tbi, by:0)

    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

    emit:
    vcfs = haplotypecaller_vcfs
    
    versions = ch_versions
}