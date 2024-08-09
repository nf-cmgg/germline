import java.nio.file.Path

// A class that contains some variables that need to be globally accesible

class GlobalVariables {
    // The available callers
    public static List availableCallers = ["haplotypecaller", "vardict"]

    public static List gvcfCallers = ["haplotypecaller"]

    public static Map<String,Path> pedFiles = [:]

}
