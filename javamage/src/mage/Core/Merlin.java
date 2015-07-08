package mage.Core;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import mage.Editor.GenbankWriter;
import mage.Editor.PlotData;
import mage.Tools.Constants;
import mage.Tools.FASTA;
import mage.optMage.Comparator;
import utils.TextFile;

/**
 *
 * Merlin is a tool for optimizing Multiplexed Automated Genome Engineering
 * oligo selection.
 *
 * The tool is designed taking into consideration Harris Wang's optMage tool as
 * well as the oligo pool considerations.
 *
 * @author Samir Ahmed
 *
 */
public class Merlin {

    // Boolean flag for enabling and disabling the compare and plot commands
    public static boolean compare = true;
    public static boolean plot = true;

    public ArrayList<Oligo> pool;

    public String targetfile = ""; //target file
    public String parameterfile = ""; //parameter file

    final static private String configFileName = "config.txt";

    /**
     * Constructor for Merlin:
     *
     * Returns an Merlin object that can be run
     *
     * @param directory	the working directory, which must have all the following
     * files
     * @param targetsFileName	the targets file, which must follow the optMage
     * input target format
     * @param parametersFileName	the parameters file, which must follow the
     * optMage input param format
     * @param genomeFileName	A fasta or ffn file
     * @throws Exception	In the event a file cannot be loaded/ parsed correctly,
     * an exception will be thrown
     */
    public Merlin(String directory, String targetsFileName, String parametersFileName, String genomeFileName) throws Exception {
        this(directory, targetsFileName, parametersFileName, genomeFileName, false);
    }

    /**
     * Constructor for Merlin:
     *
     * Returns an Merlin object that can be run
     *
     * @param directory	the working directory, which must have all the following
     * files
     * @param targetsFileName	the targets file, which must follow the optMage
     * input target format
     * @param parametersFileName	the parameters file, which must follow the
     * optMage input param format
     * @throws Exception	In the event a file cannot be loaded/ parsed correctly,
     * an exception will be thrown
     */
    public Merlin(String directory, String targetsFileName, String parametersFileName) throws Exception {
        this(directory, targetsFileName, parametersFileName, "genome.ffn", false);
    }

    //something upstream of this throws the indexoutofboundsexception.
    public void optimize() throws Exception {

        /* 
         *  Set the flags
         */
        // Set the plot condition in switches flags
        mage.Switches.Flags.plot = Merlin.plot;

        // Set the compare condition in switches flags
        mage.Switches.Flags.compare = Merlin.compare;

        /* *********************************************
         * Stage 1: Load genome /Parameters/ Switches  *
         * *********************************************/
        // First we read in the genome
        String genome = FASTA.readFFN(Oligo.Directory, Oligo.Genome);

        // Then we read and parse the parameters
        loadParameters(parameterfile);

		// The switches are public and can be turned on and off at any time
        // So there is no need to load a file
        /* ***********************************************
         * Stage 2: Load targets and Populate Oligo Pool *
         * ***********************************************/
        // Create a list of targets from by parsing the targets file
        System.out.println("Loading target file: " + targetfile);
        List<Target> targets = Target.loadTarget(targetfile);

        // Add each target in the pool to factory and push that into the pool
        for (Target tt : targets) {
            pool.add(Oligo.OligoFactory(genome, tt));
        }

        /* ****************************
         * Stage 3: Run the Heuristic *
         * ****************************/
		//
        // Optimize the pool
        mage.Core.Optimize.optimize(pool);

        /* *************************************
         * Stage 4: Return results *
         * *************************************/
    }

    /**
     * Hide all output to the error stream This will also hide any runtime
     * exceptions.
     *
     * @param Boolean enable: if true, verbose mode is enabled
     */
    public void verbose(boolean enable) {
        if (!enable) {
            System.setErr(new PrintStream(new PipedOutputStream()));
        } else {
            System.setErr(System.err);
        }
    }

    /**
     * Will search Merlin's working environment for optMage's output Oligo file
     * Will then find the corresponding oligos in the oligo folder, and save the
     * corresponding shift position
     *
     * @param filename the filename under which the desired optMage file is
     * positioned
     * @throws Exception Throws an exception if there is a problem with the IO,
     * for instance no file is found, or data in unparseable
     */
    public void compareToOptMage(String filename) throws Exception {

        // Check if the filename is valid
        if (isValidFilePath(filename)) {
            Comparator.compare(Constants.workingdirectory + filename, this.pool);
        } else {
            // Throw fnf exception
            throw new FileNotFoundException("No file called " + filename + " found");
        }

    }

    /**
     * Returns a String in Genbank format for displaying in Mage-Editor The file
     * will feature markings and the span of the oligo
     *
     * @return Genbank string
     */
    public List<String> generateGenbank() {
        ArrayList<String> list = new ArrayList<String>(pool.size());

        for (Oligo ol : pool) {
            //System.err.println("[DEBUG] Generating genbank for " + ol.name);
            GenbankWriter gw = new GenbankWriter(ol);

            list.add(gw.toString());

        }
        return list;
    }

    /**
     * Returns a List of Objects of type PlotData, which contain
     *
     * @return A List of PlotData Objects
     */
    public List<PlotData> generatePlotData() {

        ArrayList<PlotData> list = new ArrayList<PlotData>(pool.size());

        for (Oligo ol : pool) {
            PlotData pd = new PlotData(ol);
            list.add(pd);
        }

        return list;
    }

    /**
     *
     * Constructor for Merlin:
     *
     * Returns an Merlin object that can be run
     *
     * See public constructors for arguments
     *
     * @param directory
     * @param targetsFileName
     * @param parametersFileName
     * @param genomeFileName
     * @param defaultSwitches	If the default switches boolean is set to true,
     * Merlin will not search for a switches file
     * @throws Exception
     */
    private Merlin(String directory,
            String targetsFileName,
            String parametersFileName,
            String genomeFileName,
            boolean defaultSwitches) throws Exception {

        // Load the working directory
        Constants.workingdirectory = directory;

        if (!Constants.workingdirectory.endsWith("/")) {
            Constants.workingdirectory += "/";
        }
        Oligo.Directory = Constants.workingdirectory;

        // Test if the genome file exists
        if (isValidFilePath(genomeFileName)) {
            Oligo.Genome = genomeFileName;
        }

        // Test if the parameters file exists	
        if (isValidFilePath(targetsFileName)) {
            //Constants.targets = Constants.workingdirectory+ targetsFileName; 
            targetfile = directory + targetsFileName;
            Constants.targets = targetfile;
        }

        // Test if the targets file exists
        if (isValidFilePath(parametersFileName)) {
            //Constants.parameters = Constants.workingdirectory + parametersFileName; 
            parameterfile = directory + parametersFileName;
            Constants.parameters = parameterfile;
        }

        // Load the configuration files
        loadConfig(Constants.workingdirectory);

        // Create a Pool of oligos
        this.pool = new ArrayList<Oligo>();
        Oligo.resetCount();

		// Turn verbose mode off.
        //this.verbose(false);
    }

    /**
     * Helper function to load configuration text file. This file show be of the
     * following format
     * <p>
     * ____________________________ </p>
     * <p>
     * blastn < binary_path > </p>
     * <p>
     * makeblastdb  < binary_path ></p>
     * <p>
     * MFOLD < binary_path > </p>
     * <p>
     * python < binary_path ></p>
     *
     * @param workingdirectory	Constants.workingdirectory
     * @throws IOException	If the file cannot be read and parsed correctly
     */
    private void loadConfig(String workingdirectory) throws IOException {
        if (isValidFilePath(Merlin.configFileName)) {
            String configs[] = TextFile.read(workingdirectory + Merlin.configFileName).split("\n");
            for (String line : configs) {
                String[] cfg = line.split("\\s+");
                switch (cfg[0].toUpperCase()) {
                    case "BLASTN":
                        Constants.blastn = cfg[1];
                        break;
                    case "MAKEBLASTDB":
                        Constants.makeblastdb = cfg[1];
                        break;
                    case "MFOLD":
                        Constants.MFOLD = cfg[1];
                        break;
                    case "PYTHON":
                        Constants.python = cfg[1];
                        break;
                }
            }
        }
    }

    /**
     * Returns a list of with all the names of the oligos
     *
     * @return List of String
     */
    public List<String> generateNames() {
        ArrayList<String> list = new ArrayList<String>(this.pool.size());
        for (Oligo ol : pool) {
            list.add(ol.name);
        }
        return list;
    }

    /**
     * Not all the parameters in optMage are relevant here, Therefore the
     * default is given
     *
     * Parameter file requires a header (which is ignored) Parameters should be
     * oligo size dG threshold mloc_dft: (not used) "distance of mismatch to the
     * 3' end of the oligo" mloc_max: "max amount of basepair shift in the
     * oligo" cmod: number of terminal 5' phosphorothioate bonds calcreplicore:
     * automatically calculate replicore info (0=no,1=yes)
     *
     * @param filepath
     * @throws IOException
     */
    private void loadParameters(String filepath) throws IOException {

        // Read file and extract relevant data
        String[] params = TextFile.read(filepath).split("\n");
        params = params[1].split("\\s+");

        Oligo.ideal_length = Integer.parseInt(params[0]);
        Oligo.buffer_3prime = Integer.parseInt(params[3]);
        Oligo.buffer_5prime = Integer.parseInt(params[3]);

		//params[2] = Oligo.Oligo.buffer_5prime; 
    }

    /**
     * Quick Helper function for validating a file name using the
     * constants.workingdirectory
     *
     * @param filename
     * @return
     * @throws IOException
     */
    private boolean isValidFilePath(String filename) throws IOException {
        File fp = new File(Constants.workingdirectory + filename);
        if (!fp.exists()) {
            throw new IOException("Could not located file: " + Constants.workingdirectory + filename);
        }
        return true;
    }
}
