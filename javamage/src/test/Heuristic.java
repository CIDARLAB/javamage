package test;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Stack;

import mage.Core.Oligo;
import mage.Tools.FASTA;

public class Heuristic {
	
	
	public static void main(String[] args) throws Exception {
		
		Heuristic.verbose(false);
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.Blast.method = 2;
		
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  Heuristic.populate(pool);
		Heuristic.optimize(pool);
	}

	public static void optimize(ArrayList<Oligo> pool) throws Exception{
		
		// For each oligo, calculate the Blast Genome and Free Energy Values
		System.out.println("\n# Calculating Genome Homology and Free Energy Calculation [BG & DG]");
		for ( Oligo ol : pool) {
			
			ol.calc_bg();			// Calculate Blast Genome for all positions on margins
			ol.calc_dg();			// Calculate Free Energy score for all positions on margins
			
			System.out.print(ol.getOligoId()+" ");

		}
	
		// Blast all the oligos against each other
		System.out.println("\n\n# Calculating Oligo to Oligo Homology [BO]");
		for (int ii = 0; ii<(pool.size()-1); ii++ ) {
			
			// Create a list of query oligos
			ArrayList<Oligo> queries = new ArrayList<Oligo>(pool.size()-ii-1);
			
			// Generate the batch of queries
			for (int jj = ii+1; jj<pool.size(); jj++ ) {
				queries.add(pool.get(jj));
			}
			
			// Blast the queries against subject ii
			Oligo.BlastOligo(pool.get(ii),queries);
			
			// Print Oligo number
			System.out.print((ii+1)+" ");

		}
		
		// Heuristic Optimization
		System.out.println("\n\n# Heurisitic Approach\n##################");
		// Choose the oligo with the smallest mistarget score
		// Then sort by mistarget score and then repeat
		
		Stack<Oligo> stack = new Stack<Oligo>();
		stack.addAll(pool);
		for (Oligo ol: stack) {ol.reset();}
		
		int iteration = 1;
		
		while (stack.size() > 0) {

			System.out.println("\n# Iteration "+ (iteration++) +  "");
			
			// Re/Calculate BO for the entire stack
			for (Oligo ol: stack){
				ol.calc_bo();
				System.out.println("Oligo " + ol.getOligoId() + ": "+ol.scoreAt(ol.getGreedyChoice()).toString());
			}
			
			// Sort by whatever greedy-score
			Oligo.sort(stack);
			
			// Select the best choice and the repeat until the stack is empty
			Oligo greedyChoice = stack.pop();
			greedyChoice.select();
			
		}
		
		// Print the final configuration
		System.out.println("\n# Heuristic Choice");
		for (Oligo ol: pool) {
			System.out.println("Oligo "+ ol.getOligoId() + ": "+ ol.currentScore().toString() );
		}
		
	}
	
	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {
		
		Oligo.Genome = "genome0.ffn";
		Oligo.Directory = test.Constants.bo_testing;
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
	
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 250, 2,true, "o1") );
//		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 800, 2,true) );
//		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 700, 2,true) );
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 850, 2,true, "o2") );
		pool.add(Oligo.InsertionFactory(genome, "aattccgg", 650, 2,true, "o3") );
		
		return pool;
	}
	
	public static void verbose (boolean isVerbose) {
		
		if (!isVerbose){
			
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}
	

}
