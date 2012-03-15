package mage.Switches;


/**
 * Manages the Free Energy Scoring Methods
 * 
 * @author Samir Ahmed
 *
 */
public class FreeEnergy {

	public static int method = 1;
	
	/**
	 * Given a Free Energy value, it may be scored in the following ways
	 * <p>
	 * 		1 -  Offsetting the score by the threshold
	 * </p>
	 * <p>
	 * 		2 - Thresholding the score.  All value below threshold (good values) are set to zero. All other values are given a positive value
	 * </p>
	 * 
	 * @param dg_value
	 * @return
	 */
	public static Double score (Double dg_value){	

		Double score = 0.0;				// Arbitary Score

		switch (FreeEnergy.method){

		// Case 1 is 0 if greater than -12.5 or  the normalized value otherwise
		case 1: 

			if (dg_value > -12.5) { score = 0.0; }
			else { score = -12.5 + (-1.0*dg_value) ; }
			break;

			// Case 2 is just the normalized value
		case 2: score =  -12.5 + (-1.0*dg_value) ; break;
		default : System.err.println("[Switches] Invalid Scoring System Selected") ;  break;
		}
		return score;
	}
	
	/**
	 *  Setter for the threshold
	 * 
	 */
	public static boolean threshold(Double score) {
			return (score<=0.0);
	}
	
}
