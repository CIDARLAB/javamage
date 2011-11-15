package mage;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import tools.Constants;
import tools.FASTA;

public class Optimize {
	
	public static void main(String[] args) throws Exception {
		
		Optimize.verbose(true);
		mage.Switches.FreeEnergy.method = 2;
		mage.Switches.Blast.method = 1;
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  Optimize.populate(pool);
		
		// For each oligo, calculate the Blast Genome and Free Energy Values
		for ( Oligo ol : pool) {
			
			ol.calc_bg();			// Calculate Blast Genome for all positions on margins
			ol.calc_dg();			// Calculate Free Energy score for all positions on margins
			ol.calc_primaryScore();	// Calculate PrimaryScore for all positions on margins
			
			System.out.println( ol.getPrimaryScoreAsString()  );
			
		}
		
		// Blast all the oligos against each other
		for (int ii = 0; ii<pool.size(); ii++ ) {
			
			// Create a list of query oligos
			ArrayList<Oligo> queries = new ArrayList<Oligo>(pool.size()-ii-1);
			
			// Generate the batch of queries
			for (int jj = ii+1; jj<pool.size(); jj++ ) {
				queries.add(pool.get(jj));
			}
			
			// Blast the queries against subject ii
			Oligo.BlastOligo(pool.get(ii),queries);
			
		}
		
		// Calculate a starting position
		for (Oligo ol:pool) {
			ol.calc_primary_bo();
		}
		
		// Sort the pool by primary score
		Oligo.sort(pool);		
		
	} 
	
	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {
		
		String genome = FASTA.readFFN(Constants.blastdirectory,Oligo.Genome);
		
		pool.add(Oligo.InsertionFactory(genome, "GC", 190) );
		pool.add(Oligo.InsertionFactory(genome, "AT", 458) );
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
	
	public static void verbose (boolean isVerbose) {
		
		if (!isVerbose){
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}

}
