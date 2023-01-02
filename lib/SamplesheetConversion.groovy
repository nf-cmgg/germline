import groovyx.gpars.dataflow.DataflowBroadcast

class SamplesheetConversion {
    public static DataflowBroadcast convert(
        DataflowBroadcast samplesheet
    )



    //
    // Resolve Schema path relative to main workflow directory
    //
    private static String getSchemaPath(workflow, schema_filename='assets/schema_input.json') {
        return "${workflow.projectDir}/${schema_filename}"
    }
}
