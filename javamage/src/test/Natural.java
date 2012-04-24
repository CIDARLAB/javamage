package test;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Stack;

import utils.Plot;

import mage.Core.Oligo;
import mage.Tools.FASTA;

public class Natural {

	public static int  oligoNo = 50;
	public static ArrayList< ArrayList<Double[]> > bo_plots; 
	public static ArrayList< ArrayList<String  > > plot_names;

	public static void main(String[] args) throws Exception {


		Natural.verbose(false);
		//mage.Switches.FreeEnergy.method = 1;
		mage.Switches.Blast.method = 2;

		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  Natural.populateNatural(pool);

		// Optimize the pool
		Natural.OptimizeHeuristic(pool);

		// Plot the BO Values over time
		plotBO();
		
		for (Oligo ol: pool) {
			Natural.plotBG_DG(ol); 
		}
	}

	
	private static void plotBO () {
		for (int ii = 0; ii<bo_plots.size() ;ii++){
			Plot pl = new Plot();
			pl.addGraph(bo_plots.get(ii), plot_names.get(ii));
			pl.setToLines();
			pl.title("Oligo " + ii );
			pl.draw("Oligo_" + ii+"_BO" );
		}
	}

//	/**
//	 *  Plot BO Scores for a given iteration as a function of position
//	 * 
//	 */
//	private static void plotAllBO() {
//
//		// Create a new plot
//		utils.Plot pp = new utils.Plot();
//
//		// Add all the graphs
//		for (int ii = 0; ii<bo_plots.size() ;ii++){
//			pp.addGraph(bo_plots.get(ii), plot_names.get(ii));
//		}
//
//		// Set plotsytle to lines, set title and then draw
//		pp.setToLines();
//		pp.title("Blast Oligo Variation");
//		pp.draw();
//	}
	
	/**
	 * Plot the Blast Genome and Free Energy Scores as a function of position
	 * @param ol
	 */

	private static void plotBG_DG(Oligo ol){

		Double[] bgScores = ol.bgList().toArray(new Double[ol.bgList().size()]);
		Double[] dgScores  = ol.dgList().toArray(new Double[ol.dgList().size()]);

		Plot pl = new Plot();
		pl.addGraph(bgScores, "Blast Genome");
		pl.addGraph(dgScores, "Free Energy Scores"); 
		pl.setToLines();
		pl.title("Oligo : "+ol.getOligoId());
		pl.draw("Oligo_"+ol.getOligoId()+"_BG_DG");
	}

	/**
	 * Heurisitic for optimization
	 * 
	 * @param pool			An a list of Oligos 
	 * @throws Exception
	 */
	private static void OptimizeHeuristic(ArrayList<Oligo> pool) throws Exception {

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

			System.out.println("\n# Iteration "+ (iteration) +  "");

			// Re/Calculate BO for the entire stack
			for (Oligo ol: stack){
				ol.calc_bo();

				Natural.addPlot( ol.boList().toArray( new Double[ol.bgList().size()])
						, iteration ,ol.getOligoId() );

				//System.out.println("Oligo " + ol.getOligoId() + ":\t"+ol.scoreAt(ol.getGreedyChoice()).toString());
			}

			// Sort by whatever greedy-score
			Oligo.sort(stack);

			// Select the best choice and the repeat until the stack is empty
			Oligo greedyChoice = stack.pop();
			System.out.println("Oligo "+ greedyChoice.getOligoId() + ":\t"+greedyChoice.scoreAt(greedyChoice.getGreedyChoice()).toString() ); 
			greedyChoice.select();
			iteration++;
		}

		// Print the final configuration
		System.out.println("\n# Heuristic Choice");
		for (Oligo ol: pool) {
			ol.finalize();
			System.out.println("Oligo "+ ol.getOligoId() + ":\t"+ ol.currentScore().toString() );
		}

	} 

	/**
	 * Helper function for capturing the plotting
	 * 
	 * @param array
	 * @param iteration
	 * @param oligoID
	 */
	private static void addPlot(Double[] array, int iteration, int oligoID) {

		// Take array and create plot set and store
		String name = "Iteration "+iteration;
		int id = oligoID -1;
		bo_plots.get(id).add(array);
		plot_names.get(id).add(name);

	}
	
	/**
	 * populate the pool using pieces of the genome.  This means that are changing the genome as we go.
	 * To ensure that our indexing is correct, we do this 'in order', cut from the genome, remake genome
	 * and cut from next position in the genome
	 * 
	 * @param pool
	 * @return
	 * @throws Exception
	 */
	private static ArrayList<Oligo> populateNatural( ArrayList<Oligo> pool) throws Exception {

		// Load the original genome's properties
		Oligo.Genome = "genome.ffn";
		Oligo.Directory = test.Constants.naturalTestDirectory;
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
		StringBuilder gn = new StringBuilder();
		gn.append(genome);

		// Set start index
		int start = 0;
		int end = 0;
		int chunk_size = ((gn.length()/(Natural.oligoNo+1))*3)/4;

		// Create a list of target strings and list of positions
		ArrayList<String> targets = new ArrayList<String>();
		ArrayList<Integer> positions = new ArrayList<Integer>();


		// Looping and going forward
		for (int ii=0; ii < Natural.oligoNo;  ii++) {

			StringBuilder sb = new StringBuilder();
			// Define a region between the end of the last oligo and the next (genome_length/no_of_oligos) basepairs
			start 	= end + (int) Math.round((Math.random()*chunk_size));
			end 	= start + (int) Math.round(Math.random()*30);

			// Save the start positions and the target strings
			targets.add(gn.substring(start,end));
			positions.add(start);

			// Pluck out the target sequence from the genome
			sb.append(gn.substring(0,start));
			sb.append(gn.substring(end));

			// Clear the original genome and replace it with new one
			gn.setLength(0);
			gn.append(sb.toString());

			// Clear the stringbuilder buffer and garbage collect
			sb.setLength(0);
			Runtime.getRuntime().gc();
			//	System.out.println((Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) );
		}

		// Rename and write the new genome back to file
		genome = gn.toString();
		Oligo.Genome = "genomeNatural.ffn";
		FASTA.writeFFN(Oligo.Directory, Oligo.Genome, genome);

		// For each target generate an oligo
		for (int ii=0; ii<targets.size() ;ii++) {
			pool.add(Oligo.InsertionFactory(genome, targets.get(ii), positions.get(ii),2,true, Integer.toString(ii)) );
		}


		// Set up the plot List stuff
		Natural.bo_plots = new ArrayList< ArrayList<Double[]>>();
		Natural.plot_names = new ArrayList< ArrayList<String>>();
		
		for (Oligo ol : pool) {
			bo_plots.add(new ArrayList<Double[]>(ol.getMarginLength()));
			plot_names.add(new ArrayList<String>(ol.getMarginLength()));
		}

		return pool;

	}

	
	private static void verbose (boolean isVerbose) {

		if (!isVerbose){
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}

}
