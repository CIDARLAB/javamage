package mage.Switches;

import mage.Tools.BLAST.BlastResult;

/**
 * Calculate the Blast score SWITCH
 * 
 * @author Samir Ahmed
 *
 */
public class Blast {

	public static int method = 1;
	
	/**
	 * Scoring works in the following ways
	 * 
	 *  1- 	Raw some of bitscores on a blast result
	 *  2-  Weighted result of the blast result of the scores.
	 * 
	 * @param br	A blast result object for scoring
	 * @return		A double value that represents the score
	 */
	public static Double score(BlastResult br){
		Double score = 0.0;
		switch (Blast.method) {
			case 1: score = br.bitscore * Math.exp(-1.0 * br.evalue) ; break;
			case 2: score = br.qEnd-br.qStart + 1.0; break;
			default: System.err.println("[Switches] Invalid Blast Scoring System Selected") ; break;
		}
		return score;
	}
	
}
