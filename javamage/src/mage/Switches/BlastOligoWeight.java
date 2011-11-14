package mage.Switches;

import mage.Mistarget;
import mage.Oligo;

public class BlastOligoWeight {

	public	static int method = 1;
	
	
	/**
	 * Calculates the Blast Oligo Weighted score for an optimized oligo
	 * 
	 * @param ol 	Oligo to be weighted
	 * @return		Weighted Score of the Oligo based on BLAST OLIGO Mistargets
	 */
	public static Double score(Oligo ol) {

		Double score = 0.0;

		switch (mage.Switches.BlastOligoWeight.method) {
	
		// Just sum all the mistarget raw scores.
		case 1:	for(Mistarget mt : ol.valid_mt) { score += mt.getRawScore(); } break;
		default : break;
		}
		
		return score;
	}
	
}
