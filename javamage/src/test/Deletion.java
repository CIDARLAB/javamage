package test;

import mage.Core.Oligo;
import mage.Tools.Constants;
import mage.Tools.FASTA;

public class Deletion {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		mage.Switches.Blast.method = 2;
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.PrimaryScore.method =  2 ;

		try {
			String genome = FASTA.readFFN(Constants.naturalTestDirectory,"genome.ffn");
			Oligo ol = Oligo.DeletionFactory(genome, 205, 207);

			System.out.println(ol.getAsString());
			ol.calc_bg();
			System.out.println(ol.getBGasString());
			ol.calc_dg();
			System.out.println(ol.getDGasString());
			ol.calc_primaryScore();
			System.out.println("PRIMARY SCORES = " +  ol.getPrimaryScoreAsString());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}


