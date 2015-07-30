package test.Unit;

import java.io.IOException;
import java.util.List;
import mage.Core.Oligo;
import mage.Tools.FASTA;
import test.Constants;


public class OligoTesting {

	public static void main(String[] args)  {
		mage.Switches.Blast.method = 2;
		mage.Switches.FreeEnergy.method = 1;
		mage.Switches.PrimaryScore.method =  2 ;


		try {
			//testSpan();
			//testScores();
			//testPossibleOligos();
			testOneOffError(2,true);
			testOneOffError(1,true);
			testOneOffError(2,false);
			testOneOffError(1,false);
			//testReplichoreLogic();
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
		
		Oligo o1 = Oligo.MismatchFactory(genome, "CTA", 338, 342, 2, true, "t1");
		//Oligo o1 = Oligo.DeletionFactory(genome,10993,10997,2,"del");
		//Oligo o1 = Oligo.InsertionFactory(genome, "NNNNN", 1200, 2, false, "ins");
		
		System.out.println("Rep2 missense span: " + o1.span);
		System.out.println("min: " + o1.oligo_min);
		System.out.println("max: " + o1.oligo_max);
		System.out.println("margin: " + o1.getMargin());
		System.out.println("genome start: " + o1.getGenomeStart());
		System.out.println("genome end: " + o1.getGenomeEnd());
		System.out.println("target position: " + o1.target_position);
		//System.out.println("target position on genome: " + o1.x);
		System.out.println(o1.sequence);
		/*
		Oligo o2 = Oligo.InsertionFactory(genome, "NNNNN", 1200, 2, true, "ins");
		
		System.out.println("Rep2 Sense span: " + o2.span);
		System.out.println("min: " + o2.oligo_min);
		System.out.println("max: " + o2.oligo_max);
		System.out.println("margin: " + o2.getMargin());
		System.out.println("genome start: " + o2.getGenomeStart());
		System.out.println("genome end: " + o2.getGenomeEnd());
		System.out.println("target position: " + o2.target_position);
		//System.out.println("target position on genome: " + o2.x);
		System.out.println(o2.sequence);
		

		Oligo o3 = Oligo.InsertionFactory(genome, "NNNNN", 1200, 1, false, "ins");
		
		System.out.println("Rep1 Missnse span: " + o3.span);
		System.out.println("min: " + o3.oligo_min);
		System.out.println("max: " + o3.oligo_max);
		System.out.println("margin: " + o3.getMargin());
		System.out.println("genome start: " + o3.getGenomeStart());
		System.out.println("genome end: " + o3.getGenomeEnd());
		System.out.println("target position: " + o3.target_position);
		//System.out.println("target position on genome: " + o3.x);
		System.out.println(o3.sequence);
		
		Oligo o4 = Oligo.InsertionFactory(genome, "NNNNN", 1200, 1, true, "ins");
		
		System.out.println("Rep1 Sense span: " + o4.span);
		System.out.println("min: " + o4.oligo_min);
		System.out.println("max: " + o4.oligo_max);
		System.out.println("margin: " + o4.getMargin());
		System.out.println("genome start: " + o4.getGenomeStart());
		System.out.println("genome end: " + o4.getGenomeEnd());
		System.out.println("target position: " + o4.target_position);
		//System.out.println("target position on genome: " + o4.x);
		System.out.println(o4.sequence);
		*/
	}
	
	private static void testPossibleOligos() throws Exception{
		String genome = FASTA.readFFN(Constants.naturalTestDirectory,"genome.ffn");
		//Oligo o1 = Oligo.DeletionFactory(genome,10993,10997,2,"del");
		Oligo o1 = Oligo.InsertionFactory(genome, "TAG", 10993, 2, false, "ins");
		
		List<String> list = o1.getPossibleOligos();
		System.out.println("List length: " + list.size());
	}
	
	//edit for debugging
	private static void testOneOffError(int rep, boolean sense) throws Exception{
		//99 A's, TT, 10 G's, 100 T's
		String genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		System.out.println("Testing One Off Error. Rep:" + rep + " Sense:" + sense);
		//Oligo o1 = Oligo.MismatchFactory(genome, "CCC", 101, 103, 2, true, "mismatch");
		Oligo o1 = Oligo.InsertionFactory(genome, "NNN", 101, rep, sense, "ins");
		//Oligo o1 = Oligo.DeletionFactory(genome, 101, 103, 2, "del");
		
		/*System.out.println(o1.name + " span: " + o1.span);
		System.out.println("min: " + o1.oligo_min);
		System.out.println("max: " + o1.oligo_max);
		System.out.println("margin: " + o1.getMargin());
		System.out.println("genome start: " + o1.getGenomeStart());
		System.out.println("genome end: " + o1.getGenomeEnd());
		System.out.println("target position: " + o1.target_position);
		*///System.out.println("target position on genome: " + o1.x);
		System.out.println("Insert to AAATNNNTGGG");
                System.out.println(o1.sequence);
		
		//String preSequence =  genome.substring(o1.getGenomeStart()-1, 99);
		//System.out.println(preSequence);
                
                //replace the T with a C
                System.out.println("Mismatch to AAANTGGG");
                Oligo o2 = Oligo.MismatchFactory(genome, "N", 100, 101, rep, sense, "mismatch");
                System.out.println(o2.sequence);
                
                //delete the T
                System.out.println("Delete to AAATGGG");
                Oligo o3 = Oligo.DeletionFactory(genome, 100, 101, rep, "delete");
                System.out.println(o3.sequence);
                //System.out.println("deletion span: " + o3.span);
                //System.out.println("deletion genome end-start: " + (o3.getGenomeEnd()-o3.getGenomeStart()));
                //genome_end - genome_start - span  is a valid way to get deletion length

	}
	
	private static void testReplichoreLogic() throws IOException{
		String genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		int left_position = 105;
		String target = "CCGGGG";
		int replichore = 1;
		boolean sense = true;
		
		// Define the starting position on the genome)
		int genome_start = left_position + target.length() - Oligo.ideal_length + Oligo.buffer_3prime -1;
		System.out.println("Genome start: " + genome_start);
		
		// Pulls from the genome string, this start position -1 up to left position -1 +1 (not inclusive)
		String preSequence =  genome.substring(genome_start-1, left_position);
		System.out.println("Original preSeq: " + preSequence);

		// Define the ending position on the gneome
		int genome_end  = left_position + Oligo.ideal_length - Oligo.buffer_5prime -target.length();
		System.out.println("Genome end: " + genome_end);
		
		// Extract a Subsequence from the target Position to the end. THIS IS ALSO FOR A STRING i.e INDEXED FROM ZERO
		String postSequence = genome.substring(left_position , genome_end);
		System.out.println("Original postSeq: " + postSequence);

		// Take the reverse compliments of genome for -1,+1 ... ignore for -2,+2
		int splitIndex = -1;
		String reverseComp = "void";
		if (replichore == 1){

			// Get RC
			reverseComp  = mage.Tools.SequenceTools.ReverseCompliment(preSequence+postSequence);

			// Calculate index to split from and then reassign post and presequence
			splitIndex 		= postSequence.length();
			postSequence 		= reverseComp.substring(splitIndex);
			preSequence 		= reverseComp.substring(0,splitIndex);
		}

		if ( ((replichore==2) && !sense) || ((replichore==1) && sense) ) {
			target = mage.Tools.SequenceTools.ReverseCompliment(target);	
		}


		System.out.println("splitIndex: " + splitIndex);
		System.out.println("reverseComp: " + reverseComp);
		System.out.println("new post: " + postSequence);
		System.out.println("new pre: " + preSequence);
	}
}
