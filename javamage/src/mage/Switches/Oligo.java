package mage.Switches;

import java.util.List;

import mage.Core.OligoScore;


/** Class for managing Oligo and Oligo Scoring Switches 
 * 
 * @author mockingbird
 *
 */
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
	 * @return		Returns 0 if the score is not better, or 1 if it is
	 */
	public static int compare(OligoScore g1,OligoScore g2){


		int better = 0;

		switch (mage.Switches.Oligo.method_compare) {
		default: ;
		case 1: 
			
			// Tiered Lexicographic ordering system
			if (g2.FreeEnergy() < g1.FreeEnergy()) 
			{
				better++;
				
			}
			else if (g2.FreeEnergy().equals(g1.FreeEnergy()))
			{
				if (g2.BlastGenome() < g1.BlastGenome())
				{
					better++;
				}
				else if (g2.BlastGenome().equals(g1.BlastGenome()))
				{
					if (g2.BlastGenome() < g1.BlastGenome()) {
						better++;
					}
					else
					{
						better=0;
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
	
/**
 * Given an oligo, it finds the best sub oligo on the span
 * @param ol
 * @return
 */
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

	/**
	 * 
	 * OligoScore 1 gets the value of OligoScore 1+ OligoScore 2
	 * Oligo scores can be added in the following way
	 * 
	 * <p>
	 * 1 -Adding the scores together
	 * </p>
	 * @param os1 Oligo score 1
	 * @param os2 Oligo SCore 2
	 * 
	 */
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
