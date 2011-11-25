package mage.Switches;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Stack;



public class Oligo {

	public static int	method_compare 	= 1;
	public static int	method_greedy	= 1;
	/**
	 * Score comparing function
	 * 
	 * @param g1	Reference Global Score
	 * @param g2	New Global Score
	 * @return		Returns 0 if they are the same. 1 if g1 is beter or -1 if g2 is better
	 */
	public static int compare(OligoScore g1,OligoScore g2){


		int better = 1;

		switch (mage.Switches.Oligo.method_compare) {
		case 1: 
				// If g1 is worse in any return worse.
				if (	(g1.BlastGenome() 	> g2.BlastGenome() )	&&
						(g1.BlastOligo()  	> g2.BlastOligo() )		&&
						(g1.FreeEnergy() 	> g2.FreeEnergy() ) 		){
					better = 0;
				}
				
		default : better = 1 ;
		}		
		
		return better;
	}
	
	public static int greedyScore(mage.Oligo ol){
		
		return (PrimaryScore.getMinimum(ol.dgList(), ol.boList())+1);
	}
	
	
	
//	public static Stack<Integer> score(final ArrayList<Double> BG, final  ArrayList<Double> DG, final ArrayList<Double> BO, mage.Oligo ol) {
//		
//		// Depending on the priorities involved, we can choose the best solution in different ways
//		
//		switch (mage.Switches.Oligo.method_score) {
//		
//		case 1:	
//				// Here we score based on DG, then BG, then DG;
//				for (Integer ii : sortedBO) {}
//			// First we create index Lists for all Arrays
//	ArrayList<Integer> sortedBG = new ArrayList<Integer>(ol.getMarginLength()); 
//	ArrayList<Integer> sortedDG = new ArrayList<Integer>(ol.getMarginLength());
//	ArrayList<Integer> sortedBO = new ArrayList<Integer>(ol.getMarginLength());
//	
//	for (int ii = 0; ii < ol.getMarginLength() ; ii++) {
//		sortedBG.add(ii);
//		sortedDG.add(ii);
//		sortedBO.add(ii);
//	}
//		 
//	Collections.sort(sortedBG, new Comparator<Integer>() {
//		//Implementing stanard compare fucntion for sorting ascending
//		public int compare(final Integer ii, final Integer jj) {
//			return (int) ( (BG.get(ii))- (BG.get(jj)));
//		}
//	});
//	
//	Collections.sort(sortedDG, new Comparator<Integer>() {
//		//Implementing stanard compare fucntion for sorting ascending
//		public int compare(final Integer ii, final Integer jj) {
//			return (int) ( (DG.get(ii))- (DG.get(jj)));
//		}
//	});
//	
//	Collections.sort(sortedBO, new Comparator<Integer>() {
//		//Implementing stanard compare fucntion for sorting ascending
//		public int compare(final Integer ii, final Integer jj) {
//			return (int) ( (BO.get(ii))- (BO.get(jj)));
//		}
//	});
//		}
//		
//	}

}
