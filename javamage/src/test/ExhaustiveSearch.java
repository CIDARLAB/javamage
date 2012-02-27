package test;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Stack;

import mage.Core.Oligo;
import mage.Core.OligoScore;
import mage.Tools.FASTA;

public class ExhaustiveSearch {
	
	public static void main(String[] args) throws Exception {
		
		
		ExhaustiveSearch.verbose(false);
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.Blast.method = 2;
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  ExhaustiveSearch.populate(pool);
		
		
		System.out.println("\n# Calculating Genome Homology and Free Energy Calculation [BG & DG]");
		// For each oligo, calculate the Blast Genome and Free Energy Values
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
			
			System.out.print((ii+1)+" ");
		}
		
		// Heuristic Optimization
		System.out.println("\n\n# Heurisitic Approach\n##################");
		// Choose the oligo with the smallest mistarget score
		// Then sort by mistarget score and then repeat
		
		Stack<Oligo> stack = new Stack<Oligo>();
		stack.addAll(pool);
		
		int iteration = 1;
		
		while (stack.size() > 0) {

			System.out.println("\n# Iteration "+ (iteration++) +  "");
			
			// Re/Calculate BO for the entire stack
			for (Oligo ol: stack){
				ol.calc_bo();
				System.out.println("Oligo " + ol.getOligoId() + ":\t"+ol.scoreAt(ol.getGreedyChoice()).toString());
			}
			
			// Sort by whatever greedy-score
			Oligo.sort(stack);
			
			// Select the best choice and the repeat until the stack is empty
			Oligo greedyChoice = stack.pop();
			greedyChoice.select();
			
		}
		
		// Print the Heuristics configuration
		System.out.println("\n# Heuristic Choice");
		
		OligoScore globalscore = new OligoScore(); 
		
		for (Oligo ol: pool) {
		
			globalscore.add(ol.currentScore());
			System.out.println("Oligo "+ ol.getOligoId() + ":\t"+ ol.currentScore().toString() );
		}
		
		// Print the final score
		System.out.println("\n# Heuristic Score = "+globalscore.toString());
		
		
		// Exhaustive Backtracking Search
		
		
	} 
	
	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {
		
		Oligo.Genome = "genome.ffn";
		Oligo.Directory = test.Constants.bo_testing;
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
				
		for (int ii=0; ii < 20;  ii++) {
			int start =  (int) Math.round((Math.random()*genome.length()-200));
			int end = start + (int) Math.round(Math.random()*30);
			String target = genome.substring(start,end);
			
			int x = (int) Math.round((Math.random()*genome.length()-200));
			pool.add(Oligo.InsertionFactory(genome, target, x,2,true, Integer.toString(ii)) );
		}
		
		
		return pool;
		
	}
	
	public static void verbose (boolean isVerbose) {
		
		if (!isVerbose){
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}

}
