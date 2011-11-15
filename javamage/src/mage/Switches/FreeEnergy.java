package mage.Switches;

public class FreeEnergy {

	public static int method = 1;
	
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

	public static boolean threshold(Double score) {
			return (score<=0.0);
	}
	
}
