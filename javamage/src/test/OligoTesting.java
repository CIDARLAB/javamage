package test;

import mage.Oligo;
import mage.Switches;

public class OligoTesting {

	public static void main(String[] args)  {
		Switches.setBlastScoringMethod(2);
		Switches.setFreeEnergyScoringMethod(2);
		Switches.setPrimaryScoringMethod(2);
		try {	
			Oligo ol2 = new Oligo("CGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTAC",
					"GGGATTATTATTGGG","GATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGCGCGGTCACAACGTT",421,560);
		Oligo ol = Oligo.InsertionFactory(
				"CGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGCGCGGTCACAACGTT",
				"GGGATTATTATTGGG",60);
		System.out.println(ol.getOligo(1));
		System.out.println(ol2.getOligo(1));

		ol2.calc_bg();
		System.out.println(ol2.getBGasString());
		ol2.calc_dg();
		System.out.println(ol2.getDGasString());
		ol2.calc_primaryScore();
		System.out.println("PRIMARY SCORES = " +  ol2.getPrimaryScoreAsString());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
