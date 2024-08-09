
import nextflow.Nextflow
import java.nio.file.Path
import nextflow.script.WorkflowMetadata

class Pedigree {
    private Map<String,List<PedigreeEntry>> pedigrees = [:]
    Pedigree(String ped) {
        if (ped) {
            addPedContent(Nextflow.file(ped, checkIfExists:true))
        }
    }

    public void addPedContent(Path pedFile) {
        def List<String> exceptions = []
        def List<String> foundFamiliesThisPed = []
        def Integer lineCount = 0
        pedFile.readLines().each { line ->
            lineCount++
            if (line.startsWith('#')) { return }
            def PedigreeEntry entry
            try {
                entry = new PedigreeEntry(line)
            } catch(Exception e) {
                def Integer padding = 13
                exceptions.add("  - Line ${lineCount}:".padRight(padding) + "${e.getMessage().split(";;").join("\n".padRight(padding + 1))}")
                return
            }
            def String family = entry.familyId
            // Make sure families imported from other PEDs are overwritten
            if (!foundFamiliesThisPed.contains(family)) {
                pedigrees[family] = [entry]
                foundFamiliesThisPed.add(family)
            } else {
                pedigrees[family].add(entry)
            }
        }
        if (exceptions.size() > 0) {
            error("Found the following errors when parsing a PED file (${ped}):\n" + exceptions.join("\n"))
        }
    }

    public String getFamily(String sample) {
        def foundPedigree = pedigrees.find { family, entries ->
            return entries.individualId.contains(sample)
        }

        // Return samplename when no family has this sample
        def String returnValue = foundPedigree ? foundPedigree.key : sample
        return returnValue
    }

    public Map<String,Path> writePeds(WorkflowMetadata workflow) {
        def Map<String,Path> returnMap = [:]
        Nextflow.file("${workflow.workDir}/peds_${workflow.runName}").mkdir()
        pedigrees.each { String family, List<PedigreeEntry> entries ->
            def Path pedFile = Nextflow.file("${workflow.workDir}/peds_${workflow.runName}/${family}.ped")
            pedFile.text = entries.collect { it.toLine() }.join("")
            returnMap[family] = pedFile
        }
        return returnMap
    }
}

class PedigreeEntry {
    public String familyId
    public String individualId
    public String paternalId = "0"
    public String maternalId = "0"
    public String sex = "0"
    public Integer phenotype = 0

    PedigreeEntry(String pedLine) {
        def List<String> lineSplit = pedLine.split("\t") as List
        def List<String> exceptions = []

        if (lineSplit.size() < 6) {
            throw new Exception("Did not find enough values. At least 6 values are expected in this order: family_id, individual_id, paternal_id, maternal_id, sex, phenotype")
        }

        def idRegex = /^[a-zA-Z0-9_\.]+$/
        def String id

        // Family ID
        id = lineSplit[0]
        if (id ==~ idRegex) {
            familyId = id
        } else {
            exceptions.add("Invalid family ID (${id}). It should only contain these characters: a-z, A-Z, 0-9, _ and ." as String)
        }

        // Individual ID
        id = lineSplit[1]
        if (id ==~ idRegex) {
            individualId = id
        } else {
            exceptions.add("Invalid individual ID (${id}). It should only contain these characters: a-z, A-Z, 0-9, _ and ." as String)
        }

        def List<String> validMissingIDs = ["", "0"]

        // Paternal ID
        id = lineSplit[2]
        if (id ==~ idRegex) {
            paternalId = id
        } else if (!validMissingIDs.contains(id)) {
            exceptions.add("Invalid paternal ID (${id}). It should only contain these characters: a-z, A-Z, 0-9, _ and .; Use 0 if the paternal ID is missing" as String)
        }

        // Maternal ID
        id = lineSplit[3]
        if (id ==~ idRegex) {
            maternalId = id
        } else if (!validMissingIDs.contains(id)) {
            exceptions.add("Invalid maternal ID (${id}). It should only contain these characters: a-z, A-Z, 0-9, _ and .; Use 0 if the maternal ID is missing" as String)
        }

        // Sex
        def List<String> validSexValues = ["0", "1", "2"]
        id = lineSplit[4]
        sex = validSexValues.contains(id) ? id : "0"

        // Phenotype
        def List<String> validPhenotypeValues = ["-9", "0", "1", "2"]
        id = lineSplit[5]
        phenotype = validPhenotypeValues.contains(id) ? id as Integer : 0

        if (exceptions.size() > 0) {
            throw new Exception(exceptions.join(";;"))
        }
    }

    public String toLine() {
        return "${familyId}\t${individualId}\t${paternalId}\t${maternalId}\t${sex}\t${phenotype}\n" as String
    }
}
