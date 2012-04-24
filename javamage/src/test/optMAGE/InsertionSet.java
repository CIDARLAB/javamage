package test.optMAGE;

import mage.Core.Oligo;

public class InsertionSet {

	public static void main(String [] args) throws Exception
	{
		// Create a test environment for MERLIN
		String merlin_test_dir =  "/Users/mockingbird/dropbox/research/tests/insertion";
		String parametersFileName = "INPUTparam.txt"; 
		String targetFileName = "INPUTtarg.txt";
		String genomeFileName = "genome.fasta";
		
		// Create new instance of merlin object
		mage.Core.Merlin merlin = new mage.Core.Merlin(merlin_test_dir, targetFileName, parametersFileName, genomeFileName);
		merlin.verbose(false);
		// Enable plotting
		mage.Core.Merlin.plot = true;
		
		// Enable switches
		mage.Switches.Blast.method = 2;
		
		// Optimize oligos
		merlin.optimize();
		
		for (Oligo ol: merlin.pool)
		{
			System.out.println(ol.sequence);
		}
		
		merlin.compareToOptMage("OUToligos.txt");
		
		for ( Oligo ol: merlin.pool)
		{
			int op = ol.getOptMagePosition(); 
			int mp = ol.getOptimalPosition();
			double bg = ol.bgList().get(op) - ol.bgList().get(mp-1);
			double bo = ol.boList().get(op) - ol.boList().get(mp-1);
			double dg = ol.dgList().get(op) - ol.dgList().get(mp-1);
			System.out.println(dg+"\t"+bg+"\t"+bo);
		}
		
		
	}
}
