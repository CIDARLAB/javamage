import java.util.ArrayList;


public class Switches {

	private static int BlastScoringType = 1;
	private static int FreeEnergyScoringType = 1;
	private static int BlastGenomeTotallingType = 1;

	public static Double FreeEnergyScore (Double dg_value){	

		Double score = 0.0;				// Arbitary Score

		switch (Switches.FreeEnergyScoringType){

		// Case 1 is 0 if greater than -12.5 or  the normalized value otherwise
		case 1: 

			if (dg_value > -12.5) { score = 0.0; }
			else { score = dg_value+12.5 ; }
			break;

			// Case 2 is just the normalized value
		case 2: score =  -1.0*Math.abs(dg_value +12.5) ; break;
		default : System.err.println("[FreeEnergyScoring] Invalid Scoring System Selected") ;  break;
		}
		return score;
	}


	public static Double BlastScore(Double bitscore, Double evalue){
		Double score = 0.0;
		switch (Switches.BlastScoringType) {
		case 1: score = bitscore * Math.exp(-1.0 * evalue) ; break;
		case 2: score = bitscore; break;
		default: System.err.println("[BlastScoring] Invalid Blast Scoring System Selected") ; break;
		}
		return score;
	}

	public static void setBlastScoringMethod(int method){
		Switches.BlastScoringType = method;
	}

	public static void setFreeEnergyScoringMethod(int method){
		Switches.FreeEnergyScoringType = method;
	}


	public static Double BlastGenomeTotalling(ArrayList<Double> position) {
		Double total= 0.0;
		switch (Switches.BlastGenomeTotallingType) {
		case 1:  
			total = 0.0;
			for(Double ss : position) { total += ss;}		
		}
		return total;
	}

}