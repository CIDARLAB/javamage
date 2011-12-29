package test.Unit;

import java.util.HashMap;
import java.util.List;


import mage.Tools.BLAST;
import mage.Tools.BLAST.BlastResult;



public class BlastGenome {

	//Oligo is the test Oligo that we select
	public static String oligo = "TTCTACtatacgagtgcgacgtaatgcGTTGTTTGGTAGCACAAAAGTATTACCATGGTCCTAGgattacagattacagattacaAGCCT";

	//Set of the four genome with different amounts of fixed matches
	public static String[] genome = {"genome0.ffn","genome1.ffn","genome2.ffn","genome3.ffn" } ;


	public static void main (String[] args) {

		
		// Add oligo to a list for BLASTING
		HashMap<Integer,String> query = new HashMap<Integer,String>(1);
		query.put(1,oligo);

		// Blast versus Genome 0 -> 3
		for (String gg : genome) {
			BLAST bg = new BLAST(TestConstants.bgtesting_directory,gg);
			bg.setQuery(query);
			List<BlastResult> result =  bg.run();
			
			// Print out all the results Weighted followed by raw
			for (BlastResult br : result) {
				mage.Switches.Blast.method= 1;
				System.out.println("Weighted Scoring : " +mage.Switches.Blast.score(br) ) ;
				mage.Switches.Blast.method= 2;
				System.out.println("Raw Scoring : "+mage.Switches.Blast.score(br) ) ;
			}
		}
	}

}

