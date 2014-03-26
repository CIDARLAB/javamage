package test.Unit;
import java.io.IOException;
import java.util.ArrayList;

import test.Constants;
import mage.Core.Oligo;
import mage.Tools.FASTA;
import mage.Tools.OligoStats;



public class TestOligoStats {
	
	
	public static void main(String[] args){
		try {
			testGetDiversityTable();
			testGetARE();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void testGetARE() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		
		Oligo insert = Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 100, 2,true, "in1");
		double ival = OligoStats.getARE(insert);
		boolean ibool = (0.0360762 < ival && ival < 0.0360764); //0.15 * e^(-0.075*19) = 0.0360763...
		
		if (!ibool){
			System.out.println("Error in getARE(insertion): expected 0.0360763, got " + String.valueOf(ival));
		}
		
		//deletion
		Oligo delete = Oligo.DeletionFactory(genome, 200, 210, 1, "del1");
		double dval = OligoStats.getARE(delete);
		boolean dbool = (0.136466 < dval && dval < 0.136468);
		
		if (!dbool){
			System.out.println("Error in getARE(deletion): expected 0.136467, got " + String.valueOf(dval));
			System.out.println("Deletion oligo full sequence:");
			System.out.println(delete.sequence);
			System.out.println("Deletion oligo target:");
			System.out.println(delete.target);
			System.out.println("Deletion target length: " + String.valueOf(delete.target_length));
		}
		
		//mismatch
		Oligo mismatch = Oligo.MismatchFactory(genome, "ATCGA", 300, 305, 1, true, "mis1");
		double mval = OligoStats.getARE(mismatch);
		boolean mbool = (0.151514 < mval && mval < 0.151516);
		
		if (!mbool){
			System.out.println("Error in getARE(mismatch): expected 0.151515, got " + String.valueOf(mval));
		}
		
		if (mbool && dbool && ibool){
			System.out.println("getARE: PASS");
		}
	}
	
	public static void testGetCyclesForFrequency(){
		
	}
	
	public static void testGetAggregateAnyARE(){
		
	}
	
	public static void testGetAggregateSumARE(){
		
	}
	
	public static void testLociProbability() throws Exception{
		
		//String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		
		//System.out.println(genome.length());
		
		//ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		//pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		//pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );

		//1c1 * (1-(.9^1))^1 * (.9)^(1(1-2))
		double p1 = OligoStats.lociProbability(1, 1, .1, 1);
		assert (p1 == .1/.9): "Expcected lociProcability(1,1,.1,1) = .1111";
		
		
	}
	
	public static void testGetDiversityTable() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		String teststring = "abcde";
		
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );

		String res = OligoStats.getDiversityTable(pool, 10);
		System.out.println(res);
		
	}
}
