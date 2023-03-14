#!/usr/bin/env groovy

import groovy.cli.commons.CliBuilder
import groovy.cli.commons.OptionAccessor

def CliBuilder cli = new CliBuilder(usage: 'groovy merge_headers.groovy').tap {
    p(longOpt:'prefix', type: String, 'The prefix of the output file')
    v(longOpt:'vcf', type: File, 'The original VCF')
    e(longOpt:'ped_vcf', type: File, 'The VCF containing the PED header')
}

def OptionAccessor options = cli.parse(args)
def String prefix = options.p
def File vcf = options.v
def File ped_vcf = options.e

ArrayList header = []
String columns_vcf = ""
String columns_ped = ""

if(vcf.size() == 0 || ped_vcf.size() == 0){
    throw new Exception("No VCF contents detected when merging the VCF and PED headers for ${vcf} and/or ${ped_vcf}. Please contact the pipeline developer to fix this.")
}
else {
    File headerFile = new File("${prefix}.header.txt")

    vcf.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        while( line = buffered.readLine()) {
            if(line ==~ "^##.*"){
                headerFile.append(line + "\n")
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
                    headerFile.append(line + "\n")
                }
            }
            else if(line ==~ "^#CHROM.*"){
                assert line == columns_vcf : "The columns in the genotyped VCF and the VCF containing PED headers are different. (${vcf} and ${ped_vcf})"
                headerFile.append(line + "\n")
            }
        }
    }
}




