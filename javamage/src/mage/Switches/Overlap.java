package mage.Switches;

import mage.Core.Mistarget;

/** 
 * Overlap Switch 
 * @author Samir Ahmed
 *
 */
public class Overlap {

	public static int method = 1;
	
	/**
	 * Calculates the overlap between two elements
	 * 
	 * @param mt	A Mistarget object
	 * @return Vale that represents the overlaps
	 */
	public static Double score(Mistarget mt){
		
		Double score = 0.0;
		
		switch (Overlap.method) {
		case 1:
			
			// Calculate the score by scaling the overlap by the total length
			Double weighting = ( (double) (mt.overlap() ) ) / ( (double) (mt.getSequence().length()) );
			score = mt.getRawScore() * weighting;
			
		default: break;
		}
		
		return score;
		
	}
	
}
