package mage.Switches;

import java.util.ArrayList;

public class BlastGenomeTotal {

	
	public static int method = 1;

	
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
