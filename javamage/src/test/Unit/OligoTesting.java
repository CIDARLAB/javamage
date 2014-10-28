package test.Unit;

import java.io.IOException;

import test.Constants;
import mage.Core.Oligo;
import mage.Tools.FASTA;


public class OligoTesting {

	public static void main(String[] args)  {
		mage.Switches.Blast.method = 2;
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.PrimaryScore.method =  2 ;


		try {
			testSpan();
			//testScores;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void testScores(){
		try {	
			//deprecated constructors
			/*Oligo ol2 = new Oligo("CGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTAC",
					"GGGATTATTATTGGG","GATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGCGCGGTCACAACGTT",421,560,"ol2");
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
			System.out.println("PRIMARY SCORES = " +  ol2.getPrimaryScoreAsString());*/
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void testSpan() throws Exception{
		String genome = FASTA.readFFN(Constants.naturalTestDirectory,"genome.ffn");
		
		Oligo o1 = Oligo.MismatchFactory(genome, "TAG", 10993, 10997, 2, true, "t1");
		
		System.out.println("o1 span: " + o1.span);
		System.out.println("min: " + o1.oligo_min);
		System.out.println("max: " + o1.oligo_max);
	}
}
