package mage;

import java.io.IOException;
import java.util.ArrayList;

import tools.Constants;
import tools.FASTA;

public class Optimize {


	
	public static void main() throws Exception {
		
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  Optimize.populate(pool);
		
		// For each oligo, calculate the Blast Genome and Free Energy Values
		for ( Oligo ol : pool) {
			
			ol.calc_bg();
			ol.calc_dg();
			ol.calc_primaryScore();
			
		}
		
		
	} 
	
	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {
		
		String genome = FASTA.readFFN(Constants.blastdirectory,"genome.ffn");
		
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
		
		return pool;
	}

}
