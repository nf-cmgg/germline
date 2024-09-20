include { samplesheetToList } from 'plugin/nf-schema'

workflow WATCHPATH_HANDLING {

    take:
    input_samplesheet
    watchdir
    samplesheet_schema
    pedFile

    main:

    // A map with entries for every row that expects to be watched
    def watch_lines = [:]

    // Initialize a filename for the DONE file
    def watch_id = UUID.randomUUID().toString()
    def done_file = file("${watchdir}/DONE-${watch_id}")

    // A list containing all expected files
    def expected_files = [done_file.name]

    // Initialize samplesheet to keep the linter happy
    def ch_samplesheet_watched = Channel.empty()

    // Pedigree handling
    def pedigree = new Pedigree(pedFile)
    GlobalVariables.pedFiles = pedigree.writePeds(workflow)

    def errors = []

    def families = [:]
    def sample_counts = [:]
    def sample_metas = [:]

    // Determine which files to watch for
    def samplesheet_list = samplesheetToList(input_samplesheet, samplesheet_schema)
        // Do some calculations and manipulations here
        .collect { row ->
            // Watchpath logic
            def is_watch = false
            row = row.collect { input ->
                def input_name = input instanceof Path ? input.name : input as String
                if (input_name.startsWith("watch:")) {
                    is_watch = true
                    expected_files.add(input_name.replace("watch:", ""))
                    return input_name
                }
                return input
            }
            if (is_watch) {
                watch_lines[row[0].id] = row
                if (!watchdir) {
                    errors.add("Found a `watch:` prefix in the samplesheet for '${row[0].id}', but no watch directory has been set.")
                }
            }

            // Pipeline logic
            def ped = row[6]
            if (ped) {
                pedigree.addPedContent(ped)
            }
            def family = row[0].family ?: pedigree.getFamily(row[0].sample)
            def sample_id = row[0].id

            if (!families.containsKey(family)) {
                families[family] = [sample_id]
            } else if(!families[family].contains(sample_id)) {
                families[family].add(sample_id)
            }


            if (!sample_counts.containsKey(sample_id)) {
                sample_counts[sample_id] = 1
                sample_metas[sample_id] = row[0]
            } else {
                sample_counts[sample_id] += 1
                if(sample_metas[sample_id] != row[0]) {
                    def other_meta = sample_metas[sample_id]
                    def diff_keys = []
                    other_meta.each { k,v ->
                        if (v != row[0][k]) {
                            diff_keys.add(k)
                        }
                    }
                    errors.add("Found multiple entries for sample '${sample_id}' in the samplesheet with differing meta values (`${diff_keys.join(' ')}`).")
                }
            }

            def new_meta = row[0] + [family:family]
            row[0] = new_meta
            row.remove(6)
            return row
        }

    // Stop the pipeline if extra validation errors have been detected
    if (errors.size() > 0) {
        error(errors.join("\n"))
    }

    Channel.fromList(samplesheet_list).set { ch_samplesheet_all }

    if (watchdir && expected_files.size() > 1) {

        def watchdir_path = file(watchdir)
        if (!watchdir_path.exists()) {
            // Create the watchdir if it doesn't exist (this is a failsafe)
            watchdir_path.mkdir()
        }

        // Watch the watchdir for the expected files
        Channel.watchPath("${watchdir}**{${expected_files.join(',')}}", "create,modify")
            .until { file ->
                def file_name = file.name
                if (file_name == done_file.name) {
                    // Delete the done file when it's been detected and stop watching
                    done_file.delete()
                    return true
                }
                return false
            }
            .map { file ->
                // Try to find a matching ID in watch_lines
                def id = find_id(file.name, watch_lines)
                if (id == "") {
                    error("Could not find id for file '${file.name}' in the samplesheet.")
                }
                // Replace the matching watch entry with the file
                watch_lines[id] = watch_lines[id].collect { line_entry ->
                    line_entry == "watch:${file.name}" as String ? file : line_entry
                }
                return watch_lines[id]
            }
            .filter { entry ->
                def found_all_files = false
                if (!entry.any { elem -> elem.toString().startsWith("watch:") }) {
                    // Remove the entry from watch_files when all files for the current entry have been found
                    watch_lines.remove(entry[0].id)
                    found_all_files = true
                }
                if (watch_lines.size() == 0) {
                    // Create the DONE file when all files have been found
                    done_file.text = ""
                }
                // Pass through all entries where all files have been found
                return found_all_files
            }
            // Mix with all entries that didn't contain any watched files
            .mix(ch_samplesheet_all.filter { entry -> !entry.any { elem -> elem.toString().startsWith("watch:") }})
            .set { ch_samplesheet_watched }

    } else {
        ch_samplesheet_watched = ch_samplesheet_all
    }

    ch_samplesheet_watched
        .map { row ->
            row[0] = row[0] + [
                    family_samples:families[row[0].family].sort(false).join(","),
                    duplicate_count:sample_counts[row[0].id]
                ]
            return row
        }
        .set { samplesheet }

    emit:
    samplesheet

}

// Find the ID of a file in a map with sample IDs as keys
def find_id(file_name, file_map) {
    def lastDotIndex = file_name.lastIndexOf(".")
    if (lastDotIndex == -1) {
        return ""
    }
    def id = file_name.substring(0, lastDotIndex)
    if (file_map.containsKey(id)) {
        return id
    }
    return find_id(id, file_map)
}
