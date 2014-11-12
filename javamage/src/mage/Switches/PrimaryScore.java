package mage.Switches;

import java.util.ArrayList;
import java.util.List;

public class PrimaryScore {

	// Scoring Method
	public	static int	method= 2;
	
	/**
	 * Function that determines how the primary score is calculated
	 * 
	 * 
	 * @param 	bg_scores	An list of Blast Genome Score, generated after calling calc_bg
	 * @param 	dg_scores	A list of Free Energy Scores, generated after calling calc_dg 
	 * @return	an integer index of the start position of the optimal oligo based on the primary score
	 * 
	 */
	//TODO: Find if this is ever invoked outside of test classes
	public static int score(ArrayList<Double> bg_scores, ArrayList<Double> dg_scores) {

		int primary_position = -1;

		switch (PrimaryScore.method){

		case 1:  primary_position = getMinimum(bg_scores,dg_scores); break;
		case 2:  primary_position = getMinimum(dg_scores,bg_scores); break;
		}

		return primary_position+1;
	}

	/**
	 * Function that determines how the primary score is calculated
	 * 
	 * 
	 * @param 	bg_scores	An list of Blast Genome Score, generated after calling calc_bg
	 * @param 	dg_scores	A list of Free Energy Scores, generated after calling calc_dg 
	 * @param	bo_scores	A list of Blast scores against the oligo pool
	 * @return	an integer index of the start position of the optimal oligo based on the primary score
	 * 
	 */
	//TODO: Confirm you can do this- are BO scores created after this is set?
	public static int score(ArrayList<Double> bg_scores, ArrayList<Double> dg_scores, ArrayList<Double> bo_scores) {

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
		//find the lowest score in the main list and its first position
		for (Double dd : list) {
			if (dd < min) { index = counter; min = dd;}
			counter++;
		}

		ArrayList<Double> possible = new ArrayList<Double>();		
		ArrayList<Integer> indicies = new ArrayList<Integer>();

		counter =0;		
		//get the scores in list 2 at the positions that match the lowest
		//score of list 1
		for (Double dd : list) {

			if (dd.equals(min) ) {
				possible.add(list2.get(counter)); 
				indicies.add(counter);  
			}
			counter ++;
		}

		Double min2 = possible.get(0);
		counter = 0;
		//find the lowest score/position in list2's allowed positions
		for (Double dd : possible) {
			if ( dd < min2 ) { index = indicies.get(counter); min2 = dd;}
			counter++;
		}


		return index;
	}
	
	/**
	 * Function for determing the minimum values of list 2 given a list 1
	 * 
	 * Remember that is indexed from Zero
	 * 
	 * 
	 * @param list	List with more priority
	 * @param list2 List with second priority
	 * @param list3 list with third priority
	 * @return		Index number of the minimum value
	 * 
	 */
	public static int getMinimum(List<Double> list, List<Double> list2, List<Double> list3){
		Double min = list.get(0);
		Integer counter = 0;
		Integer index = 0;
		
		//find the lowest score in the main list and its first position
		for (Double dd : list) {
			if (dd < min) { index = counter; min = dd;}
			counter++;
		}

		ArrayList<Double> possible = new ArrayList<Double>();		
		ArrayList<Integer> indicies = new ArrayList<Integer>();

		counter =0;
		//get the scores in list 2 at the positions that match the lowest
		//score of list 1
		for (Double dd : list) {

			if (dd.equals(min) ) {
				possible.add(list2.get(counter)); 
				indicies.add(counter);  
			}
			counter ++;
		}

		Double min2 = possible.get(0);
		counter = 0;
		
		//get the scores in list 3 at the positions that match the lowest
		//score of list 1/2
		ArrayList<Double> possible2 = new ArrayList<Double>();		
		ArrayList<Integer> indicies2 = new ArrayList<Integer>();
		counter = 0;
		for (Double dd: possible){
			if (dd.equals(min2)){
				possible2.add(list3.get(counter));
				indicies2.add(counter);
			}
			counter++;
		}

		//find the lowest score/position in list3's allowed positions
		Double min3 = possible2.get(0);
		counter = 0;
		for (Double dd: possible2){
			if (dd < min3) {index = indicies2.get(counter); min3 = dd;}
			counter++;
		}
		
		return index;
	}
	
}
