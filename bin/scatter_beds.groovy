#!/usr/bin/env groovy

import groovy.cli.commons.CliBuilder
import groovy.cli.commons.OptionAccessor

def CliBuilder cli = new CliBuilder(usage: 'groovy scatter_beds.groovy').tap {
    p(longOpt:'prefix', type: String, 'The prefix of the output files')
    b(longOpt:'bed', type: File, 'The original BED')
    s(longOpt:'size', type: Integer, 'The minimal size the BEDs should be')
}

def OptionAccessor options = cli.parse(args)
def String prefix = options.p
def File bed = options.b
def Integer scatterSize = options.s

def String fileName = bed.toString().toLowerCase()

if (fileName.endsWith("bed")) {
    bed.withReader {
        Integer fileCount = 1
        String status = "new"
        Integer currentSize = 0
        File outputFile
        while( line = it.readLine() ) {
            if( line ==~ /^#.*$/) {continue} // skip comments
            if(status == "new"){
                outputFile = new File("${prefix}_${fileCount}.bed")
                status = "busy"
            }
            outputFile.append("${line}\n")
            def ArrayList entry = line.tokenize("\t")
            currentSize = currentSize + (entry[2].toInteger() - entry[1].toInteger() + 1)
            if( currentSize >= scatterSize ){
                currentSize = 0
                fileCount++
                status = "new"
            }
        }
    }
}
else {
    throw new Exception("${bed} has an unknown file extension. It has to be a BED file.")
}

