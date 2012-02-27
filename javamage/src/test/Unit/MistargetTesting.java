package test.Unit;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import mage.Core.Mistarget;
import mage.Core.Oligo;
import test.Constants;
import mage.Tools.FASTA;

public class MistargetTesting {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {

		MistargetTesting.verbose(false);
		long start = System.currentTimeMillis();
		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  MistargetTesting.populate(pool);


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
		
		System.out.println("\nMistarget Count "+Mistarget.mistarget_collection.size());
		
		for (Oligo ol: pool) { 
			ol.calc_bg();
			ol.calc_dg();
			ol.calc_primaryScore();
			ol.set(ol.getPrimaryPosition());
			
			// Now calculate the weighted_bo score of every oligo (Done by a switch)
			ol.calc_primary_bo();
			
			System.out.println("Oligo " + ol.getOligoId() + " | Primary Postion: " + ol.getPrimaryPosition() + " | Weighted BO Score:" + ol.getWeightedBOScore());
		}
		
		// Sort the oligos based on their wieghted scores
		Oligo.sort(pool);
		
		
		for (Oligo ol : pool) {
		System.out.println("ID: "+ol.getOligoId()+ "	Score: "+ol.getWeightedBOScore());		
		ol.calc_bo();	
		}

		System.out.println("\nRuntime: "+(System.currentTimeMillis()-start)+" ms");

	} 

	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {

		//Oligo.Genome = "genome2.ffn";
		String genome = FASTA.readFFN(Constants.blastdirectory,Oligo.Genome);
		
		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );

//		pool.add(Oligo.InsertionFactory(genome, "atgtgtaagtacgtgcaccagt", 250));
//		pool.add(Oligo.InsertionFactory(genome, "atgatgcgcgcgatagcgagatcgac", 800));
		
		return pool;
	}

	public static void verbose (boolean isVerbose) {

		if (!isVerbose){
			System.setErr( new PrintStream( new PipedOutputStream() ) );
		}
	}


}

