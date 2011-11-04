package mage;
import java.util.ArrayList;
import java.util.List;

import tools.BLAST.BlastResult;


public class Switches {

	private static int SumMistargetMethod = 1;
	private static int PrimaryScoringMethod = 1;
	private static int BlastScoringType = 1;
	private static int FreeEnergyScoringType = 1;
	private static int BlastGenomeTotallingType = 1;

	public static Double FreeEnergyScore (Double dg_value){	

		Double score = 0.0;				// Arbitary Score

		switch (Switches.FreeEnergyScoringType){

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


	public static Double BlastScore(BlastResult br){
		Double score = 0.0;
		switch (Switches.BlastScoringType) {
		case 1: score = br.bitscore * Math.exp(-1.0 * br.evalue) ; break;
		case 2: score = br.qEnd-br.qStart + 1.0; break;
		default: System.err.println("[Switches] Invalid Blast Scoring System Selected") ; break;
		}
		return score;
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

	
	public static int PrimaryScore(ArrayList<Double> bg_scores, ArrayList<Double> dg_scores) {
		
		int primary_position = -1;
		
		switch (Switches.PrimaryScoringMethod){
		
		case 1:  primary_position = getMinimum(bg_scores,dg_scores); break;
		case 2:  primary_position = getMinimum(dg_scores,bg_scores); break;
		}
		
		return primary_position;
	}

	public static int getMinimum(List<Double> list, List<Double> list2){
		Double min = list.get(0);
		Integer counter = 0;
		Integer index = 0;
		for (Double dd : list) {
			if (dd < min) { index = counter; min = dd;}
			counter++;
		}
		
		ArrayList<Double> possible = new ArrayList<Double>();		
		ArrayList<Integer> indicies = new ArrayList<Integer>();
		
		counter =0;
		for (Double dd : list) {

			if (dd.equals(min) ) {
				possible.add(list2.get(counter)); 
				indicies.add(counter);  
			}
			counter ++;
		}
		
		Double min2 = possible.get(0);
		counter = 0;
		for (Double dd : possible) {
			if ( dd < min2 ) { index = indicies.get(counter); min2 = dd;}
			counter++;
		}
		
		
		return index;
	}
	
	public static Double SumMistargetScores(ArrayList<Double> scores) {
		
		double total=0.0;
		switch (Switches.SumMistargetMethod) {
		
		case 1: for(Double ss: scores) {total+=ss;} break;
			
		}
		return total;
	}
	
	public static void setBlastScoringMethod(int method){
		Switches.BlastScoringType = method;
	}

	public static void setFreeEnergyScoringMethod(int method){
		Switches.FreeEnergyScoringType = method;
	}
	
	public static void setPrimaryScoringMethod(int method){
		Switches.PrimaryScoringMethod = method;
	}

}


