package test.Unit;

import java.util.ArrayList;


import mage.Core.Oligo;
import test.Constants;
import mage.Tools.FASTA;
import utils.TextFile;

public class testDGScoring {

	
	private static final String directory = "/Users/mockingbird/dropbox/research/optimization/testFiles/testDGScoreMethod/";

	/**
	 * Short Static Function for printing and plotting results of DG calculations
	 * 
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
	

		String genome = FASTA.readFFN(Constants.blastdirectory,Oligo.Genome);
		
		System.out.println(genome.length());
		
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );


		Integer count = 1; 
		mage.Switches.FreeEnergy.method = 2;
		// Print files to the specified directory with the output results
		for ( Oligo ol: pool ){
			ol.calc_dg();
			String out = "# Oligo Number "+count+" - Target Sequence: "+ol.getTarget()+"\n"+ol.getDGasString(); 
			String filename = "Oligo_"+count.toString() + ".txt";
			TextFile.write(directory,filename, out);
			System.out.println("[TEST] Generated Outfile " + count.toString());
			count++;
		}
		

	}

}
