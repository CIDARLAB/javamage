package mage.Switches;

import mage.Core.Mistarget;
import mage.Core.Oligo;

public class BlastOligo {

	public	static int method = 1;


	/**
	 * Calculates the Blast Oligo Weighted score for an optimized oligo
	 * 
	 *  Scoring works in the following way/s 
	 *  <p>
	 *  1 - Sum the mistarget scores  
	 *  </p>
	 * @param ol 	Oligo to be weighted
	 * @return		Weighted Score of the Oligo based on BLAST OLIGO Mistargets
	 */
	public static Double score(Oligo ol) {

		Double score = 0.0;

		switch (mage.Switches.BlastOligo.method) {

		// Just sum all the mistarget raw scores.
		case 1:	for(Mistarget mt : ol.valid_mt) {
			score += /*mt.getRawScore()*/mt.score(); 
			}
		break;
		default : break;
		}

		return score;
	}

}
