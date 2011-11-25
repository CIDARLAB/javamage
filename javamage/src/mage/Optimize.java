package mage;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Stack;

import tools.Constants;
import tools.FASTA;

public class Optimize {
	
	public static void main(String[] args) throws Exception {
		
		Optimize.verbose(true);
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.Blast.method = 2;
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  Optimize.populate(pool);
		
		// For each oligo, calculate the Blast Genome and Free Energy Values
		for ( Oligo ol : pool) {
			
			ol.calc_bg();			// Calculate Blast Genome for all positions on margins
			ol.calc_dg();			// Calculate Free Energy score for all positions on margins
		}
		
		// Blast all the oligos against each other
		for (int ii = 0; ii<(pool.size()-1); ii++ ) {
			
			// Create a list of query oligos
			ArrayList<Oligo> queries = new ArrayList<Oligo>(pool.size()-ii-1);
			
			// Generate the batch of queries
			for (int jj = ii+1; jj<pool.size(); jj++ ) {
				queries.add(pool.get(jj));
			}
			
			// Blast the queries against subject ii
			Oligo.BlastOligo(pool.get(ii),queries);
			
		}
		
		// Heuristic Optimization
		// Choose the oligo with the smallest mistarget score
		// Then sort by mistarget score and then repeat
		
		Stack<Oligo> stack = new Stack<Oligo>();
		stack.addAll(pool);
		
		while (stack.size() > 0) {

			// Re/Calculate BO for the entire stack
			for (Oligo ol: stack){
				ol.calc_bo();
				ol.reset();
				System.out.println(ol.getBOasString());	
			}
			
			// Sort by whatever greedy-score
			Oligo.sort(stack);
			
			// Select the best choice and the repeat until the stack is empty
			Oligo greedyChoice = stack.pop();
			greedyChoice.select();
			
		}
			
		
	} 
	
	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {
		
		Oligo.Genome = "genome0.ffn";
		Oligo.Directory = Constants.bo_testing;
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
	
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 250) );
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 850) );
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 650) );
		
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
