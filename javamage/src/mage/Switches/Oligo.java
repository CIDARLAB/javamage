package mage.Switches;

import mage.Core.OligoScore;

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
	
	public static int greedyScore(mage.Core.Oligo ol){
		
		return (PrimaryScore.getMinimum(ol.dgList(), ol.boList())+1);
	}
	

}
