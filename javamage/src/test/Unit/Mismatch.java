package test.Unit;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import test.Constants;

public class Mismatch {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		mage.Switches.Blast.method = 2;
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.PrimaryScore.method =  2 ;

		try {
			String genome = FASTA.readFFN(Constants.naturalTestDirectory,"genome.ffn");
			
			//int random_position = (int) Math.floor(Math.random()*100000 + 100);
			Oligo ol = mage.Core.Oligo.MismatchFactory(
							genome, 
							"ggaattaaccaa",
							//genome.substring( random_position, random_position + 30 ), 
							205, 
							207, 
							2,
							true,
							"olTest");
			
			
			System.out.println(ol.getAsString());
			ol.calc_bg();
			System.out.println(ol.getBGasString());
			ol.calc_dg();
			System.out.println(ol.getDGasString());

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
