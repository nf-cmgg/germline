import groovyx.gpars.dataflow.DataflowBroadcast
import groovy.json.JsonSlurper
import org.yaml.snakeyaml.Yaml
import java.nio.file.Path
import nextflow.Nextflow
import nextflow.Channel


class SamplesheetConversion {
    public static DataflowBroadcast convert(
        Path samplesheetFile,
        Path schemaFile
    ) {

        def Map schema = (Map) new JsonSlurper().parseText(schemaFile.text)
        def Map schemaFields = schema.get("properties")
        def ArrayList allFields = schemaFields.keySet().collect()
        def ArrayList requiredFields = schema.get("required")

        def String fileType = getFileType(samplesheetFile)
        def String delimiter = fileType == "csv" ? "," : fileType == "tsv" ? "\t" : null
        def DataflowBroadcast samplesheet

        if(fileType == "yaml"){
            samplesheet = Channel.fromList(new Yaml().load((samplesheetFile.text)))
        }
        else {
            samplesheet = Channel.value(samplesheetFile).splitCsv(header:true, strip:true, sep:delimiter)
        }

        // Field checks + returning the channels
        def Map uniques = [:]
        def Boolean headerCheck = true
        def Integer sampleCount = 0

        return samplesheet.map({ row ->

            sampleCount++

            // Check the header once for CSV/TSV and for every sample for YAML
            if(headerCheck) {
                def ArrayList rowKeys = row.keySet().collect()
                def ArrayList differences = allFields.plus(rowKeys)
                differences.removeAll(allFields.intersect(rowKeys))

                def String yamlInfo = fileType == "yaml" ? " for sample ${sampleCount}." : ""

                def ArrayList samplesheetDifferences = rowKeys.intersect(differences)
                if(samplesheetDifferences.size > 0) {
                    throw new Exception("[Samplesheet Error] The samplesheet contains following unwanted field(s): ${samplesheetDifferences}${yamlInfo}")
                }

                def ArrayList requiredDifferences = requiredFields.intersect(differences)
                if(requiredDifferences.size > 0) {
                    throw new Exception("[Samplesheet Error] The samplesheet requires '${requiredFields.join(",")}' as header field(s), but is missing these: ${requiredDifferences}${yamlInfo}")
                }

                if(fileType in ["csv", "tsv"]) {
                    headerCheck = false
                }
            }

            // Check required dependencies
            def Map dependencies = schema.get("dependentRequired")
            if(dependencies) {
                for( dependency in dependencies ){
                    if(row[dependency.key] != "" && row[dependency.key]) {
                        def ArrayList missingValues = []
                        for( value in dependency.value ){
                            if(row[value] == "" || !row[value]) {
                                missingValues.add(value)
                            }
                        }
                        if (missingValues) {
                            throw new Exception("[Samplesheet Error] ${dependency.value} field(s) should be defined when '${dependency.key}' is specified, but  the field(s) ${missingValues} are/is not defined.")
                        }
                    }
                }
            }

            def Map meta = [:]
            def ArrayList output = []

            for( field : schemaFields ){
                def String key = field.key
                def String regexPattern = field.value.pattern && field.value.pattern != '' ? field.value.pattern : '^.*$'
                def String metaNames = field.value.meta

                def String input = row[key]

                if((input == null || input == "") && key in requiredFields){
                    throw new Exception("[Samplesheet Error] Sample ${sampleCount} does not contain an input for required field '${key}'.")
                }
                else if(field.value.unique){
                    if(!(key in uniques)){uniques[key] = []}
                    if(input in uniques[key] && input){
                        throw new Exception("[Samplesheet Error] The '${key}' value needs to be unique. '${input}' was found twice in the samplesheet.")
                    }
                    uniques[key].add(input)
                }
                else if(!(input ==~ regexPattern) && input != '' && input) {
                    throw new Exception("[Samplesheet Error] The '${key}' value for sample ${sampleCount} does not match the pattern '${regexPattern}'.")
                }

                if(metaNames) {
                    for(name : metaNames.tokenize(',')) {
                        meta[name] = (input != '' && input) ? input.replace(' ', '_') : field.value.default ?: null
                    }
                }
                else {
                    def inputFile = (input != '' && input) ? Nextflow.file(input) : field.value.default ? Nextflow.file(field.value.default) : []
                    if( inputFile != [] && !inputFile.exists() ){
                        throw new Exception("[Samplesheet Error] The '${key}' file (${input}) for sample ${sampleCount} does not exist.")
                    }
                    output.add(inputFile)

                }
            }
            output.add(0, meta)
            return output
        })

    }

    private static String getFileType(
        Path samplesheetFile
    ) {
        def String extension = samplesheetFile.getExtension()
        if (extension in ["csv", "tsv", "yml", "yaml"]) {
            return extension == "yml" ? "yaml" : extension
        }

        def String header = getHeader(samplesheetFile)

        def Integer commaCount = header.count(",")
        def Integer tabCount = header.count("\t")

        if ( commaCount == tabCount ){
            throw new Exception("[Samplesheet Error] Could not derive file type from ${samplesheetFile}. Please specify the file extension (CSV, TSV, YML and YAML are supported).")
        }
        if ( commaCount > tabCount ){
            return "csv"
        }
        else {
            return "tsv"
        }
    }

    private static String getHeader(
        Path samplesheetFile
    ) {
        def String header
        samplesheetFile.withReader { header = it.readLine() }
        return header
    }
}
