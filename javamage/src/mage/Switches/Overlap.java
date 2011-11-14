package mage.Switches;

import mage.Mistarget;

public class Overlap {

	public static int method = 1;
	
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
