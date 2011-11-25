package mage.Switches;

import java.util.ArrayList;
import java.util.List;

public class PrimaryScore {

	// Scoring Method
	public	static int	method= 1;
	
	/**
	 * Function that determines how the primary score is calculated
	 * 
	 * 
	 * @param 	bg_scores	An list of Blast Genome Score, generated after calling calc_bg
	 * @param 	dg_scores	A list of Free Energy Scores, generated after calling calc_dg 
	 * @return	an integer index of the start position of the optimal oligo based on the primary score
	 * 
	 */
	public static int score(ArrayList<Double> bg_scores, ArrayList<Double> dg_scores) {

		int primary_position = -1;

		switch (PrimaryScore.method){

		case 1:  primary_position = getMinimum(bg_scores,dg_scores); break;
		case 2:  primary_position = getMinimum(dg_scores,bg_scores); break;
		}

		return primary_position+1;
	}

	/**
	 * Function for determing the minimum values of list 2 given a list 1
	 * 
	 * Remember that is indexed from Zero
	 * 
	 * 
	 * @param list	List with more priority
	 * @param list2 List with second priority
	 * @return		Index number of the minimum value of L1 and L2
	 * 
	 */
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
	
}
