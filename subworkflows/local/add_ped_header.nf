//
// ADD PED HEADER
//

include { RTGTOOLS_PEDFILTER                  } from '../../modules/nf-core/rtgtools/pedfilter/main'
include { BCFTOOLS_REHEADER                   } from '../../modules/nf-core/bcftools/reheader/main'

workflow ADD_PED_HEADER {
    take:
        vcfs                 // channel: [mandatory] [ meta, vcfs ] => The post-processed VCFs
        somalier_samples_tsv // channel: [mandatory] [ meta, samples_tsv ] => The samples TSV retrieved from SOMALIER_RELATE

    main:

    ch_versions         = Channel.empty()

    //
    // Remove extra columns from the samples TSV and convert to a VCF header
    //

    somalier_samples_tsv
        .map { meta, samples_tsv ->
            convert_to_ped(samples_tsv)
        }

    RTGTOOLS_PEDFILTER(
        somalier_samples_tsv
    )

    ch_versions = ch_versions.mix(RTGTOOLS_PEDFILTER.out.versions)

    //
    // Create the new headers
    //

    vcfs
        .join(RTGTOOLS_PEDFILTER.out.output)
        .multiMap { meta, vcf, ped_vcf ->
            header: [ meta, add_ped_to_header(vcf, ped_vcf) ]
            files: [ meta.id, meta, vcf]
        }
        .set { headers }

    headers.header
        .collectFile()
        { meta, header ->
            [ "${meta.id}.header.txt", header.join("\n") ]
        }
        .map { file ->
            id = file.baseName.tokenize(".")[0..-1].join(".")
            [ id, file ]
        }
        .join(headers.files)
        .map { id, header, meta, vcf ->
            [ meta, vcf, header ]
        }
        .dump(tag:'bcftools_reheader_input', pretty:true)
        .set { bcftools_reheader_input }

    BCFTOOLS_REHEADER(
        bcftools_reheader_input,
        []
    )

    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    emit:
    ped_vcfs = BCFTOOLS_REHEADER.out.vcf
    versions = ch_versions
}

def convert_to_ped(samples_tsv) {

    ArrayList new_lines = []

    for(line : samples_tsv.readLines()){
        new_line = line.tokenize("\t")[0..5].join("\t")
        new_lines.add(new_line)
    }

    samples_tsv.text = new_lines.join("\n")
}

def add_ped_to_header(vcf, ped_vcf) {

    ArrayList header = []
    String columns_vcf = ""
    String columns_ped = ""

    if(vcf.size() == 0 || ped_vcf.size() == 0){
        println("WARNING: No VCF contents detected when merging the VCF and PED headers for ${vcf} and/or ${ped_vcf} - ignore this when using stub runs")
        header.add("##STUB")
    }
    else {
        vcf.withInputStream {
            InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            BufferedReader buffered = new BufferedReader(decoder)
            while( line = buffered.readLine()) {
                if(line ==~ "^##.*"){
                    header.add(line)
                }
                else if(line ==~ "^#CHROM.*"){
                    columns_vcf = line
                }
                else {
                    break
                }
            }
        }

        ped_vcf.withInputStream {
            InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            BufferedReader buffered = new BufferedReader(decoder)
            while( line = buffered.readLine()) {
                if(line ==~ "^##.*"){
                    if (!(line ==~ "^##file.*")){
                        header.add(line)
                    }
                }
                else if(line ==~ "^#CHROM.*"){
                    columns_ped = line
                }
            }
        }
    }



    assert columns_ped == columns_vcf : "The columns in the genotyped VCF and the VCF containing PED headers are different. (${vcf} and ${ped_vcf})"
    header.add(columns_vcf)
    return header

}
