package test.Unit;

import mage.Core.Oligo;

public class OligoTesting {

	public static void main(String[] args)  {
		mage.Switches.Blast.method = 2;
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.PrimaryScore.method =  2 ;
		try {	
			Oligo ol2 = new Oligo("CGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTAC",
					"GGGATTATTATTGGG","GATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGCGCGGTCACAACGTT",421,560,"ol2");
//			Oligo ol = Oligo.InsertionFactory(
//					"CGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGCGCGGTCACAACGTT",
//					"GGGATTATTATTGGG",60);
//			System.out.println(ol.getOligo(1));
//			System.out.println(ol2.getOligo(1));

			ol2.calc_bg();
			System.out.println(ol2.getBGasString());
			ol2.calc_dg();
			System.out.println(ol2.getDGasString());
			//ol2.calc_primaryScore();
		//	System.out.println("PRIMARY SCORES = " +  ol2.getPrimaryScoreAsString());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
