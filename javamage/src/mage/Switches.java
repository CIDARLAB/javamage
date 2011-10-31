package mage;
import java.util.ArrayList;
import java.util.List;


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
	
	public static void setPrimaryScoringMethod(int method){
		Switches.PrimaryScoringMethod = method;
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

	private static int getMinimum(List<Double> list, List<Double> list2){
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
	
	
	public static void main( String[] args) {
		ArrayList<Double> l1= new ArrayList<Double>();
		ArrayList<Double> l2= new ArrayList<Double>();
		
		l1.add(94.0);
		l1.add(97.0);
		l1.add(94.0);
		l1.add(97.0);
		l1.add(98.0);
		l1.add(14.0);
		l1.add(14.0);
		l1.add(34.0);

		l2.add(0.0);
		l2.add(0.0);
		l2.add(0.0);
		l2.add(0.0);
		l2.add(0.0);
		l2.add(4.0);
		l2.add(1.0);
		l2.add(0.0);
		
		int index = getMinimum(l1,l2);
		System.out.println( index +  " index gives "+ l1.get(index) ) ;
		
	}


	public static Double SumMistargetScores(ArrayList<Double> scores) {
		
		double total=0.0;
		switch (Switches.SumMistargetMethod) {
		
		case 1: for(Double ss: scores) {total+=ss;} break;
			
		}
		return total;
		
	}
}


