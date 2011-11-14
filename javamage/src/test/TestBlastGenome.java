package test;

import java.util.ArrayList;

import mage.Oligo;
import tools.Constants;
import tools.FASTA;
import utils.TextFile;


public class TestBlastGenome {

	private final static String directory = "/Users/mockingbird/dropbox/research/optimization/testFiles/testBlastGenomeScores/";
	
	
	/**
	 *  Test function for BLAST Genome Scores
	 *  
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {

		
		
		String genome = FASTA.readFFN(Constants.blastdirectory,"genome.ffn");
		
		System.out.println(genome.length());
		
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190) );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458) );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408) );
		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349) );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599) );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658) );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900) );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029) );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426) );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592) );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925) );

		Integer count = 1; 
		mage.Switches.Blast.method = 2;
		// Print files to the specified directory with the output results
		for ( Oligo ol: pool ){
			ol.calc_bg();
			String out = "# Oligo Number "+count+" - Target Sequence: "+ol.getTarget()+"\n"+ol.getBGasString(); 
			String filename = "Oligo_"+count.toString() + ".txt";
			TextFile.write(directory,filename, out);
			System.out.println("[TEST] Generated Outfile " + count.toString());
			count++;
		}

	}
}


