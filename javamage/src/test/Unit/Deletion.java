package test.Unit;

import java.io.IOException;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import test.Constants;

public class Deletion {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		mage.Switches.Blast.method = 2;
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.PrimaryScore.method =  2 ;

		try {
			//testScores();
			testPosition();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void testScores() throws Exception{
		String genome = FASTA.readFFN(Constants.naturalTestDirectory,"genome.ffn");
		Oligo ol = Oligo.DeletionFactory(genome, 205, 207, 2,"deletionTest");

		System.out.println(ol.getAsString());
		ol.calc_bg();
		System.out.println(ol.getBGasString());
		ol.calc_dg();
		System.out.println(ol.getDGasString());
		ol.calc_primaryScore();
		System.out.println("PRIMARY SCORES = " +  ol.getPrimaryScoreAsString());
	}
	
	//ideal behavior: left_postion is the position (starting the count at 1) of the first deleted base. Right_position is the position of the first base kept after the deletion
	//right_position - left_position = deletion length
	private static void testPosition() throws Exception{
		String genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGGGGNCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
		//that's 100 A's, TGGGGN, 100 C's
		//Now, delete the Gs
		
		Oligo rep1 = Oligo.DeletionFactory(genome, 102, 105, 1, "Rep1");
		System.out.println(rep1.sequence);
		Oligo rep2 = Oligo.DeletionFactory(genome, 102, 105, 2, "Rep2");
		System.out.println(rep2.sequence);
	}
}


