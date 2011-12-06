package mage.Switches;

import java.util.List;

import mage.Core.OligoScore;

public class Oligo {

	public static int	method_compare 	= 1;
	public static int	method_greedy	= 1;
	public static int	method_add		= 1;

	// Constant weighting
	public static double	dg_constant	= 100.0;
	public static double	bg_constant = 1.0;
	public static double	bo_constant = 1.0;
	/**
	 * Score comparing function
	 * 
	 * @param g1	Reference Global Score
	 * @param g2	New Global Score
	 * @return		Returns 0 if they are the same. 1 if g1 is beter or -1 if g2 is better
	 */
	public static int compare(OligoScore g1,OligoScore g2){


		int better = 0;

		switch (mage.Switches.Oligo.method_compare) {
		default: ;
		case 1: 
			// Lexicographical ordering
			if (g2.FreeEnergy() <= g1.FreeEnergy()) {
				if (g2.BlastGenome() <= g1.BlastGenome() ){
					if (g2.BlastOligo() <= g1.BlastOligo()) {
						better= 1;
					}
				}
			}
			break;
		case 2: 
			// Weighted summing
			double d1 = dg_constant*g1.FreeEnergy() + bg_constant*g1.BlastGenome() + bo_constant*g1.BlastOligo();
			double d2 = dg_constant*g2.FreeEnergy() + bg_constant*g2.BlastGenome() + bo_constant*g2.BlastOligo();
			if (d2 > d1) {
				better =1;
			}
			break;


		}		

		return better;
	}

	public static int greedyScore(mage.Core.Oligo ol){

		List<OligoScore> list = ol.getScores();
		int choice=-1;

		switch (method_greedy) {

		default: ;
		case 1 :
			//  Compare scores and return the best

			OligoScore best =list.get(0);
			int counter = 1;
			choice = counter;

			// Loop through and get the minimum 
			for (OligoScore os : list) {
				if ( os.isBetterThan(best) ) {
					best = os;
					choice = counter;
				}
				counter++;
			}
			break;
		}

		return choice;
	}

	public static void add(OligoScore os1, OligoScore os2) {

		switch (Oligo.method_add){

		default:
		case 1: 	os1.BlastGenome(os1.BlastGenome() 	+ os2.BlastGenome() );
					os1.BlastOligo(	os1.BlastOligo() 	+ os2.BlastOligo()	);
					os1.FreeEnergy(	os1.FreeEnergy()	+ os2.FreeEnergy()	);
		break;
		}
	}


}
