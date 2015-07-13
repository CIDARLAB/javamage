package mage.Tools;

public final class Constants {

    // Defaults
    public static String makeblastdb = "c:/Program Files/NCBI/blast-2.2.28+/bin/makeblastdb";
    public static String blastn = "c:/Program Files/NCBI/blast-2.2.28+/bin/blastn";
    public static String MFOLD = "/usr/local/bin/mfold";
    //point to the Python executable, not Biopython
    public static String python = "python";

    public static String workingdirectory = System.getProperty("user.dir") + System.getProperty("file.separator");

    public static String targets = "";
    public static String parameters = "";
    public static String switches = "";

}
