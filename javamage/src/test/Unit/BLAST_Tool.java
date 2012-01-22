package test.Unit;

import java.util.HashMap;

import mage.Tools.BLAST;
import test.Constants;

public class BLAST_Tool {

	public static void main(String[] args) {

		String directory =  Constants.blastdirectory;
		String subjectFFN = "genome.ffn";

		HashMap <Integer,String> seqMap = new HashMap<Integer,String>();
		seqMap.put(12,"CAAATGTGCAGCATACGTATTTGCTCGGCGTGCTTGGTCTCTCGTACTTCTCCTGGAGATCAAGGAAATGTTTCTTGTCCAAGCGGACAG");
		seqMap.put(144,"AGCTCTCGATAGACGACTACGGGAAAATACGATCAGGACTCGGGACTACGATATACGACAGAAAGACGTACGATTTACGACTGCGCGCGA");

		BLAST bg = new BLAST(directory,subjectFFN);
		bg.setQuery(seqMap);
		bg.run();

	} 
	
	
}
