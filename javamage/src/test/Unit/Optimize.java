package test.Unit;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import mage.Core.Oligo;
import mage.Tools.FASTA;

public class Optimize {
	
	
	public static int  oligoNo = 10;
	
	/**
	 * Test function for optimize oligo
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {


		verbose(false);
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.Blast.method = 2;

		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  populateNatural(pool);

		// Optimize the pool
		mage.Core.Optimize.optimize(pool);
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
		int chunk_size = ((gn.length()/(oligoNo+1))*3)/4;

		// Create a list of target strings and list of positions
		ArrayList<String> targets = new ArrayList<String>();
		ArrayList<Integer> positions = new ArrayList<Integer>();


		// Looping and going forward
		for (int ii=0; ii < oligoNo;  ii++) {

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

		return pool;

	}
	
	
	private static void verbose (boolean isVerbose) {

		if (!isVerbose){
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}
	
}
