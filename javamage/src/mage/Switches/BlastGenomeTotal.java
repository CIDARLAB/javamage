package mage.Switches;

import java.util.ArrayList;


/** 
 * Class for managing how the blast genome score is totaled
 * 
 * @author Samir Ahmed
 *
 */
public class BlastGenomeTotal {

	
	public static int method = 1;

	/**
	 * scoring works in the following way/s
	 * 
	 * 	1- Sum of all the doubles in the list
	 * 
	 * @param position An array list of the scores. 
	 * @return The total blast genome score
	 */
	public static Double score(ArrayList<Double> position) {
		Double total= 0.0;
		switch (BlastGenomeTotal.method){
		case 1:  
			total = 0.0;
			for(Double ss : position) { total += ss;}		
			default : break;
		}
		return total;
	}

}
