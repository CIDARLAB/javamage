
public class Switches {

	private static int BlastScoring = 1;
	private static int FreeEnergyScoring = 1;

	public static Double FreeEnergyScore (Double dg_value){	

		Double score = 0.0;				// Arbitary Score

		switch (Switches.FreeEnergyScoring){

		// Case 1 is 0 if greater than -12.5 or  the normalized value otherwise
		case 1: 

			if (dg_value > -12.5) { score = 0.0; }
			else { score = dg_value-12.5 ; }
			break;

			// Case 2 is just the normalized value
		case 2: score =  dg_value +12.5 ; break;
		default : System.err.println("[FreeEnergyScoring] Invalid Scoring System Selected") ;  break;
		}
		return score;
	}


	public static Double BlastScore(Double bitscore, Double evalue){
		Double score = 0.0;
		switch (Switches.BlastScoring) {
		case 1: score = bitscore * Math.exp(-1.0 * evalue) ; break;
		case 2: score = bitscore; break;
		default: System.err.println("[BlastScoring] Invalid Blast Scoring System Selected") ; break;
		}
		return score;
	}

	public static void setBlastScoringMethod(int method){
		Switches.BlastScoring = method;
	}

	public static void setFreeEnergyScoringMethod(int method){
		Switches.FreeEnergyScoring = method;
	}

}