import java.nio.file.Path

// A class that contains some variables that need to be globally accesible

class GlobalVariables {
    // The available callers
    public static List availableCallers = ["haplotypecaller", "vardict", "elprep"]

    public static List gvcfCallers = ["haplotypecaller", "elprep"]

    public static List bamCallers = ["elprep", "vardict"]

    public static Map<String,Path> pedFiles = [:]

}
