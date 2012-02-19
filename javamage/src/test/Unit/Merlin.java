package test.Unit;

public class Merlin {

	/**
	 * Test Function
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		// Create a test environment for MERLIN
		String merlin_test_dir =  "/Users/mockingbird/dropbox/research/optimization/merlin";
		String parametersFileName = "INPUTparam.txt"; 
		String targetFileName = "INPUTtarg.txt";
		String genomeFileName = "genome.fasta";
		
		// Create new instance of merlin object
		mage.Core.Merlin merlin = new mage.Core.Merlin(merlin_test_dir, targetFileName, parametersFileName, genomeFileName);
		
		// Enable plotting
		mage.Core.Merlin.plot = true;
		
		// Enable switches
		mage.Switches.Blast.method = 2;
		
		// Optimize oligos
		merlin.optimize();
		
		merlin.compareToOptMage("OUToligos.txt");
		
	}

}
