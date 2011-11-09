package test;

import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import mage.Mistarget;
import mage.Oligo;
import tools.BLAST;
import tools.BLAST.BlastResult;
import tools.Constants;
import tools.FASTA;

public class MistargetTesting {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {


		// Create a collection of oligos and populate it
		ArrayList<Oligo> pool = new ArrayList<Oligo> ();
		pool  =  MistargetTesting.populate(pool);
		
		Oligo.BlastOligo(pool.get(0), pool.get(1));
		
//		BLAST bg = new BLAST("/Users/mockingbird/dropbox/research/optimization/mt_testing/","genome.ffn");
//
//		HashMap<Integer,String> query = new HashMap<Integer,String>(1);
//		int count = 0;
//
//		for (Oligo ol: pool){
//			System.out.println("Testing Oligo: " +ol.getSequenceAsString().toString());
//			query.put(count++, ol.getSequenceAsString().toString());
//		}
//		bg.setQuery(query);
//		List<BlastResult> result =  bg.run();
//		System.err.println(result.size());
//		
//		ArrayList<Mistarget> mt = new ArrayList<Mistarget> (result.size());
//		for (BlastResult br: result){
//			mt.add(new Mistarget(br,0,br.oligoID));
//		}
		
//		
//		for (Mistarget mm: mt) {
//			
//		}
		

	} 

	private static ArrayList<Oligo> populate( ArrayList<Oligo> pool) throws Exception {

		String genome = FASTA.readFFN(Constants.blastdirectory,Oligo.Genome);

		pool.add(Oligo.InsertionFactory(genome, "GC", 190) );
		//		pool.add(Oligo.InsertionFactory(genome, "AT", 458) );
		//		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408) );
		//		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349) );
		//pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599) );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658) );
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

