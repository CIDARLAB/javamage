package test.Unit;

import java.util.List;

public class GenbankTester {
	
	public static void main(String [] args) throws Exception {
		
//		String genome = FASTA.readFFN(Constants.blastdirectory,Oligo.Genome);
//		Oligo ol = Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "testOligo");
//		GenbankWriter  test= new GenbankWriter(ol);
//		System.out.println(test.toString() );
//		
//		System.out.println(test.genbankSequenceFormat(ol.getSequenceAsString()));
		
		
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
		
		try{
			List<String> list = merlin.generateGenbank();
		
			for ( String ss : list )
			{
				System.out.println("\n----------\n"+ss);
			}
		}
		catch(Exception e){System.err.println("Error encountered generating Genbank file");}
		
		
	}
}
