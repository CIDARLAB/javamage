package test.Unit;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import test.Constants;
import mage.Core.Oligo;
import mage.Tools.FASTA;
import mage.Tools.OligoStats;



public class TestOligoStats {
	
	public static void main(String[] args){
		try {
			//testGetDiversityTable2();
			//testGetARE();
			//testLociProbability();
			//testGetAggregateAnyARE();
			//testGetAggregateSumARE();
			testOligoLength();
			//checkOligoStructure();
			//testGetARE2();
			//testGetCumulativeDiversityTable();
			//testGetDiscreteDiversityTable();
			//testProbabilitiesForAllCycles();
			//testProbabilitiesAfterCycles();
		} catch (Exception e) {
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
		Oligo delete = Oligo.DeletionFactory(genome, 200, 212, 1, "del1");
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
			//System.out.println("Spans (insert, delete, mismatch");
			//System.out.println(String.valueOf(insert.span));
			//System.out.println(String.valueOf(delete.span));
			//System.out.println(String.valueOf(mismatch.span));
		}
	}
	
	public static void testGetARE2() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();
		
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1785136, 2, false, "pps"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 3079657, 2, false, "tktA"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 2739172, 2, false, "aroF"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 784855, 1, true, "aroG"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1786458, 2, true, "aroH"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1772709, 2, true, "aroD"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 3428860, 2, false, "aroE"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1771812, 2, true, "ydiB"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 3517086, 2, false, "aroK"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 405628, 1, true, "aroL"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 958034, 1, true, "aroA"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 2445495, 2, false, "aroC"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1320970, 1, false, "trpE"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1319408, 1, false, "trpD"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1317809, 1, false, "trpC"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1315246, 1, false, "trpA"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 1316439, 1, false, "trpB"));
		pool.add(Oligo.InsertionFactory(genome, "TAG", 3484141, 2, true, "crp"));
		
		for (Oligo oligo : pool){
			System.out.println(oligo.name + ": " + OligoStats.getARE(oligo));
		}
		System.out.println("getAggregateAnyARE: " + OligoStats.getAggregateAnyARE(pool));
		System.out.println("getAggregateSumARE: " + OligoStats.getAggregateSumARE(pool));
	}
	
	public static void testGetCyclesForFrequency(){
		
	}
	
	public static void testGetAggregateAnyARE() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		
		Oligo insert = Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 100, 2,true, "in1");
		Oligo delete = Oligo.DeletionFactory(genome, 200, 212, 1, "del1");
		Oligo mismatch = Oligo.MismatchFactory(genome, "ATCGA", 300, 305, 1, true, "mis1");
		
		ArrayList<Oligo> oligos = new ArrayList<Oligo>();
		oligos.add(insert);
		oligos.add(delete);
		oligos.add(mismatch);
		oligos.add(insert);
		oligos.add(delete);
		oligos.add(mismatch);

		System.out.println("Aggregate Any ARE: " + String.valueOf(OligoStats.getAggregateAnyARE(oligos)));
	}
	
	public static void testGetAggregateSumARE() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		
		Oligo insert = Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 100, 2,true, "in1");
		Oligo delete = Oligo.DeletionFactory(genome, 200, 212, 1, "del1");
		Oligo mismatch = Oligo.MismatchFactory(genome, "ATCGA", 300, 305, 1, true, "mis1");
		
		ArrayList<Oligo> oligos = new ArrayList<Oligo>();
		oligos.add(insert);
		oligos.add(delete);
		oligos.add(mismatch);

		System.out.println("Aggregate Sum ARE: " + String.valueOf(OligoStats.getAggregateSumARE(oligos)));	
	}
	
	public static void testLociProbability() throws Exception{

		//1c1 * (1-(.9^1))^1 * (.9)^(1(1-2))
		double p1 = OligoStats.lociProbability(2, 2, .2, 2);
		System.out.println("Expected lociProcability(2,2,0.2,2) = .1296, got " + String.valueOf(p1));
		
		double p2 = OligoStats.lociProbability(3,2,.05,2);
		assert (p2 > .35911 && p2 < .35912): "Expected lociProbability(3,2,.05,2)";
		System.out.println("Expected lociProcability(3,2,.05,2) = .0316, got " + String.valueOf(p2));
		
		double p3 = OligoStats.lociProbability(5,5,.5,5);
		System.out.println("Expected lociProcability(5,5,.5,5) = .8532, got " + String.valueOf(p3));
	}
	
	public static void testGetDiversityTable() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );//.036
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );//.045
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );//.071
		pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );//.096
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );//.071
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo"));//.018
		pool.add(Oligo.InsertionFactory(genome, "GACGATATGAAC", 7900, 2,true, "oligo") );//.066
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );//.082
		pool.add(Oligo.InsertionFactory(genome, "ATTTGACCT", 15426, 2,true, "oligo") );//.082
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );//.027
		//pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );

		String res = OligoStats.getDiversityTable(pool, 10);
		System.out.println("Diversity table 1:");
		System.out.println(res);
		System.out.println();
		
		System.out.println("Any ARE: " + OligoStats.getAggregateAnyARE(pool));
		System.out.println("Sum ARE: " + OligoStats.getAggregateSumARE(pool));
	}
	
	public static void testGetDiversityTable2() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();
		pool.add(Oligo.MismatchFactory(genome, "TAG", 5000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 15000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 25000, 5003, 2, true, "oligo"));
		/*pool.add(Oligo.MismatchFactory(genome, "TAG", 35000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 45000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 55000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 65000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 75000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 85000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 95000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 105000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 115000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 125000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 135000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 145000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 155000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 165000, 5003, 2, true, "oligo"));
		pool.add(Oligo.MismatchFactory(genome, "TAG", 175000, 5003, 2, true, "oligo"));
		*/
		String res = OligoStats.getDiversityTable(pool, 10);
		System.out.println("Diversity table 1:");
		System.out.println(res);
		System.out.println();
		
		System.out.println("Any ARE: " + OligoStats.getAggregateAnyARE(pool));
		System.out.println("Sum ARE: " + OligoStats.getAggregateSumARE(pool));

	}
	
	//confirm that Oligos are in fact 90bp long
	public static void testOligoLength() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		Oligo insert = Oligo.InsertionFactory(genome, "AAAAAAAAAAAAAA", 190, 2,true, "oligo");
		Oligo delete = Oligo.DeletionFactory(genome, 200, 212, 1, "del1");
		Oligo mismatch = Oligo.MismatchFactory(genome, "ATCGA", 300, 305, 1, true, "mis1");

		System.out.println("Insert:" + insert.sequence);
		System.out.println("Deletion:" + delete.sequence);
		System.out.println("Mismatch:" + mismatch.sequence);
		
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
		
		for (Oligo ol:pool){
			System.out.println(ol.sequence.length());
		}
		
	}
	
	//output all the information about an oligo that you can
	public static void checkOligoStructure() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		//Oligo insert = Oligo.InsertionFactory(genome, "AAAAAAAAAAAAAA", 190, 2,true, "oligo");
		Oligo insert = Oligo.DeletionFactory(genome, 90, 100, 2, "oligo");
		Map<String, Object> imap = new HashMap<String, Object>();
		imap.put("sequence", insert.sequence);
		imap.put("sequence.length", insert.sequence.length());
		imap.put("span", insert.span);
		imap.put("target", insert.target);
		imap.put("target_length", insert.target_length);
		imap.put("target_position", insert.target_position);
		//imap.put("x", insert.x);
		imap.put("valid_mt", insert.valid_mt);
		imap.put("getBioBegin()", insert.getBioBegin());
		imap.put("getBioEnd()", insert.getBioEnd());
		imap.put("getGenomeStart()", insert.getGenomeStart());
		imap.put("getGenomeEnd()", insert.getGenomeEnd());
		imap.put("getOligoMinimum()", insert.getOligoMinimum());
		imap.put("getOligoMaxmimum()", insert.getOligoMaximum());
		imap.put("getOptimized()", insert.getOptimized());
		imap.put("getOptimized().length()", insert.getOptimized().length());
		imap.put("getOptimizedStart()", insert.getOptimizedStart());
		imap.put("getOptimizedEnd()", insert.getOptimizedEnd());
		imap.put("getOptmagePosition()", insert.getOptMagePosition());
		imap.put("getPossibleOligos()", insert.getPossibleOligos());
		//imap.put("getSource()" , insert.getSource());
		imap.put("getTarget()", insert.getTarget());
		imap.put("buffer_3prime", insert.buffer_3prime);
		imap.put("buffer_5prime", insert.buffer_5prime);
		
		ArrayList<Oligo> pool = new ArrayList<Oligo>();
		pool.add(insert);
		mage.Core.Optimize.optimize(pool);  
		
		System.out.println("Insertion oligo contents:");
		List<String> keys = new ArrayList<String>(imap.keySet());
		Collections.sort(keys);
		for (String key : keys){
			System.out.println(key + " : " + imap.get(key).toString());
		}
		
	}
	
	private static void testGetCumulativeDiversityTable() throws Exception {
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		/*pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );
		*/
		
		System.out.println("Cumulative Diversity");
		System.out.println(OligoStats.getCumulativeDiversityTable(pool, 10));
	}
	
	private static void testGetDiscreteDiversityTable() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		/*pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );
		*/
		System.out.println("Discrete Diversity");
		System.out.println(OligoStats.getDiscreteDiversityTable(pool, 10));
	}
	
	private static void testProbabilitiesForAllCycles() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		/*pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );
		*/
		System.out.println("Probability For All Cycles");
		System.out.println(OligoStats.probabilitiesForAllCycles(pool, 10));
	}
	
	private static void testProbabilitiesAfterCycles() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
		/*pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACG", 7900, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "AT", 15426, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );
		pool.add(Oligo.InsertionFactory(genome, "GATTACGACCAGT", 1514925, 2,true, "oligo") );
		*/
		System.out.println("Probability After 1 Cycle");
		System.out.println(OligoStats.probabilitiesAfterCycles(pool, 1));
		System.out.println("Probability After 2 Cycles");
		System.out.println(OligoStats.probabilitiesAfterCycles(pool, 2));
	}
}
