package mage.Core;

import java.io.File;
import java.io.IOException;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import mage.Tools.Constants;
import mage.Tools.FASTA;
import mage.optMage.Comparator;
import utils.TextFile;

/**
 * 
 * Merlin is a tool for optimizing Multiplexed Automated Genome Engineering
 * oligo selection.
 * 
 * The tool is designed taking into consideration Harris Wang's optMage tool
 * as well as the oligo pool considerations.
 * 
 * @author Samir Ahmed
 *
 */
public class Merlin{
	
	// Boolean flag for enabling and disabling the compare and plot commands
	public static boolean compare = true;
	public static boolean plot = true;
	
	public ArrayList<Oligo> pool ;
	
	final static private String configFileName= "config.txt";
	
	/**
	 * Constructor for Merlin: 
	 * 
	 * Returns an Merlin object that can be run
	 * 
	 * @param directory				the working directory, which must have all the following files
	 * @param targetsFileName		the targets file, which must follow the optMage input target format
	 * @param parametersFileName	the parameters file, which must follow the optMage input param format
	 * @param genomeFileName		A fasta or ffn file
	 * @throws Exception			In the event a file cannot be loaded/ parsed correctly, an exception will be thrown
	 */
	
	public Merlin (String directory, String targetsFileName, String parametersFileName, String genomeFileName) throws Exception {
		this( directory, targetsFileName, parametersFileName, genomeFileName ,false);
	}
	
	
	/**
	 * Constructor for Merlin: 
	 * 
	 * Returns an Merlin object that can be run
	 * 
	 * @param directory				the working directory, which must have all the following files
	 * @param targetsFileName		the targets file, which must follow the optMage input target format
	 * @param parametersFileName	the parameters file, which must follow the optMage input param format
	 * @throws Exception			In the event a file cannot be loaded/ parsed correctly, an exception will be thrown
	 */
	public Merlin (String directory, String targetsFileName, String parametersFileName ) throws Exception { 
		this( directory, targetsFileName, parametersFileName, "genome.ffn" ,false);
	}
	
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
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
		
		// Then we read and parse the parameters
		loadParameters(Constants.parameters);
		
		// The switches are public and can be turned on and off at any time
		// So there is no need to load a file
		
		/* ***********************************************
		 * Stage 2: Load targets and Populate Oligo Pool *
		 * ***********************************************/
		
		// Create a list of targets from by parsing the targets file
		List<Target> targets = Target.loadTarget(Constants.targets);
		
		// Add each target in the pool to factory and push that into the pool
		for ( Target tt : targets) {
			pool.add( Oligo.OligoFactory(genome, tt) );
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
	 * Hide all output to the error stream
	 * This will also hide any runtime exceptions.
	 * 
	 * @param Boolean enable:  if true, verbose mode is enabled
	 */
	public void verbose(boolean enable){
		if (!enable){
			System.setErr( new PrintStream( new PipedOutputStream() ) );	
		}
		else{
			System.setErr( System.err);
		}
	}
	
	public void compareToOptMage( String filename) throws Exception{
		
		// Check if the filename is valid
		if ( isValidFilePath(filename) ) { 
			Comparator.compare(Constants.workingdirectory+filename, this.pool);
		}
		
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
	 * @param defaultSwitches		If the default switches boolean is set to true, Merlin will not search for a switches file
	 * @throws Exception
	 */
	private Merlin (String directory,
			String targetsFileName, 
			String parametersFileName,  
			String genomeFileName,
			boolean defaultSwitches) throws Exception {
		
		// Load the working directory
		Constants.workingdirectory = directory;
		
		if ( !Constants.workingdirectory.endsWith("/") ) { Constants.workingdirectory += "/"; }		
		Oligo.Directory = Constants.workingdirectory;
		
		// Test if the genome file exists
		if ( isValidFilePath( genomeFileName )) { Oligo.Genome = genomeFileName; }
		
		// Test if the parameters file exists	
		if ( isValidFilePath( targetsFileName)) { Constants.targets = Constants.workingdirectory+ targetsFileName; }
		
		// Test if the targets file exists
		if ( isValidFilePath( parametersFileName)) { Constants.parameters = Constants.workingdirectory + parametersFileName; }	
		
		// Load the configuration files
		loadConfig(Constants.workingdirectory);
		
		// Create a Pool of oligos
		this.pool = new ArrayList<Oligo>();
				
		// Turn verbose mode off.
		this.verbose(false);
	} 
	
	/**
	 * Helper function to load configuration text file. This file show be of the following format
	 * <p>____________________________ </p>
	 * <p>blastn < binary_path > </p>
	 * <p>makeblastdb  < binary_path ></p>
	 * <p>MFOLD < binary_path > </p>
	 * 
	 * @param workingdirectory	Constants.workingdirectory
	 * @throws IOException		If the file cannot be read and parsed correctly
	 */
	private void loadConfig(String workingdirectory) throws IOException {
		if (isValidFilePath(Merlin.configFileName)) {
			String configs[] = TextFile.read(workingdirectory + Merlin.configFileName).split("\n");
			Constants.blastn = configs[1].split("\\s+")[1];
			Constants.makeblastdb = configs[0].split("\\s+")[1];
			Constants.MFOLD = configs[2].split("\\s+")[1];
		}
	}

	/**
	 * Not all the parameters in optMage are relevant here,
	 * Therefore the default is given 
	 * @param filepath
	 * @throws IOException
	 */
	private void loadParameters(String filepath) throws IOException {
		
		// Read file and extract relevant data
		String[] params = TextFile.read(filepath).split("\n");
		params = params[1].split("\\s+");
		
		Oligo.ideal_length = Integer.parseInt(params[0]);
		Oligo.buffer_3prime = Integer.parseInt(params[1]);
		
		//params[2] = Oligo.Oligo.buffer_5prime; 
		
	}
	
	/**
	 * Quick Helper function for validating a file name using the constants.workingdirectory
	 * 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	private boolean isValidFilePath( String filename) throws IOException { 
		File fp =  new File( Constants.workingdirectory +filename ) ;
		if (! fp.exists()) { throw new IOException("Could not located file: " + Constants.workingdirectory + filename ); }
		return true;
	}
}
