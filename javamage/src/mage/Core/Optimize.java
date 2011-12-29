//package mage.Core;
//
//import java.io.IOException;
//import java.io.PipedOutputStream;
//import java.io.PrintStream;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Stack;
//
//import mage.Tools.Constants;
//
//import utils.TextFile;
//
//public class Optimize {
//
//
//	public ArrayList< ArrayList<Double[]> > bo_plots; 
//	public ArrayList<ArrayList<String   > > plot_names;
//
//
//	public static void main(String[] args) {
//		
//	}
//	
//	private static void setSwitches(){
//		
//		Constants.workingdirectory =System.getProperty("user.dir")+"/";
//		TextFile.read("working directory");
//	} 
//	
//	public Optimize (List<Oligo> pool) throws IOException {
//		this.bo_plots = new ArrayList< ArrayList<Double[]> >();
//		this.plot_names = new ArrayList< ArrayList<String> >();
//		
//		blast(pool);
//	}
//	
//	private void blast(List<Oligo> pool) throws IOException {
//
//		System.out.println("\n# Calculating Genome Homology and Free Energy Calculation [BG & DG]");
//
//		// For each oligo, calculate the Blast Genome and Free Energy Values
//		for ( Oligo ol : pool) {
//
//			ol.calc_bg();			// Calculate Blast Genome for all positions on margins
//			ol.calc_dg();			// Calculate Free Energy score for all positions on margins
//
//			System.out.print(ol.getOligoId()+" ");
//		}
//
//		// Blast all the oligos against each other
//		System.out.println("\n\n# Calculating Oligo to Oligo Homology [BO]");
//		for (int ii = 0; ii<(pool.size()-1); ii++ ) {
//
//			// Create a list of query oligos
//			ArrayList<Oligo> queries = new ArrayList<Oligo>(pool.size()-ii-1);
//
//			// Generate the batch of queries
//			for (int jj = ii+1; jj<pool.size(); jj++ ) {
//				queries.add(pool.get(jj));
//			}
//
//			// Blast the queries against subject ii
//			Oligo.BlastOligo(pool.get(ii),queries);
//
//			System.out.print((ii+1)+" ");
//		}
//
//	}
//
//	public void heuristic(List<Oligo> pool) throws Exception{
//		// Heuristic Optimization
//		System.out.println("\n\n# Heurisitic Approach\n##################");
//
//
//		// Choose the oligo with the smallest mistarget score
//		// Then sort by mistarget score and then repeat	
//		Stack<Oligo> stack = new Stack<Oligo>();
//		stack.addAll(pool);
//
//		int iteration = 1;
//
//		while (stack.size() > 0) {
//
//			System.out.println("\n# Iteration "+ (iteration) +  "");
//
//			// Re/Calculate BO for the entire stack
//			for (Oligo ol: stack){
//				ol.calc_bo();
//
//				addPlot( ol.boList().toArray( new Double[ol.bgList().size()])
//						, iteration ,ol.getOligoId() );
//
//				System.out.println("Oligo " + ol.getOligoId() + ":\t"+ol.scoreAt(ol.getGreedyChoice()).toString());
//			}
//
//			// Sort by whatever greedy-score
//			Oligo.sort(stack);
//
//			// Select the best choice and the repeat until the stack is empty
//			Oligo greedyChoice = stack.pop();
//			greedyChoice.select();
//			iteration++;
//		}
//
//		// Print the final configuration
//		System.out.println("\n# Heuristic Choice");
//		for (Oligo ol: pool) {
//
//			System.out.println("Oligo "+ ol.getOligoId() + ":\t"+ ol.currentScore().toString() );
//		}
//	}
//
//	private void addPlot(Double[] array, int iteration, int oligoID) {
//
//		// Take array and create plot set and store
//		String name = "Iteration "+iteration;
//		int id = oligoID -1;
//		this.bo_plots.get(id).add(array);
//		this.plot_names.get(id).add(name);
//
//	}
//
//
//	private static void verbose (boolean isVerbose) {
//
//		if (!isVerbose){
//			System.setErr( new PrintStream( new PipedOutputStream() ) );
//		}
//	}
//
//}
