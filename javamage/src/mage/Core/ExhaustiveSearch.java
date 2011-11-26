package mage.Core;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import mage.Tools.Constants;
import mage.Tools.FASTA;



public class ExhaustiveSearch {

	/**
	 * Conducts and Exhaustive Search, calculating every possible BO, BG, DG combination
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {

		ExhaustiveSearch.verbose(true);
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.Blast.method = 2;
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  ExhaustiveSearch.populate(pool);

		// For each oligo, calculate the Blast Genome and Free Energy Values
		for ( Oligo ol : pool) {
			ol.calc_bg();			// Calculate Blast Genome for all positions on margins
			ol.calc_dg();			// Calculate Free Energy score for all positions on margins
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
		
			// Exhaustive Search
			search(pool,0);
		
	}
		
	public static void search(List<Oligo> pool,int oligoID) {
		
	}
	

	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {

		Oligo.Genome = "genome0.ffn";
		Oligo.Directory = Constants.bo_testing;
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);

		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 250) );
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 850) );
		//		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408) );
		//		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349) );
		//		pool.add(Oligo.InsertionFactory(genome, "GC", 190) );
		//		pool.add(Oligo.InsertionFactory(genome, "AT", 458) );
		//		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408) );
		//		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349) );
		//		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599) );
		//		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658) );
		//		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900) );
		//		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029) );
		//		pool.add(Oligo.InsertionFactory(genome, "AT", 15426) );
		//		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592) );
		//		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925) );

		return pool;
	}

	public static void verbose (boolean isVerbose) {

		if (!isVerbose){
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}

}

