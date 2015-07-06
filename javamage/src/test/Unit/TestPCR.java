package test.Unit;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import test.Constants;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;
import mage.Tools.PCR;

public class TestPCR {

	//Isaacs et al supplement put the modified base at the end (3') of the forward primer.
	//Should fix forward primer design to put modified bases towards the back of the sequence, not the start
	//confirmed from same paper that the forward primer is not reversed or complementary
	
	//here's a chunk of sequence we use from the E coli genome
	//remember substrings are built with 0 at the first index
	//80 TTCTGAACTG
	//90 GTTACCTGCC
	//100 GTGAGTAAAT
	//110 TAAAATTTTA

	public static void main(String[] args) throws Exception{
		testGetUnmodifiedForwardPrimerFromOligo();
		testGetModifiedForwardPrimerFromOligo();
		//testGetMASCPCRPrimersForOligo();
		testPCRbyReplicatingResults();
                testGeneratePrimerSets();
	}
	
	public static void testGetUnmodifiedForwardPrimerFromOligo() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2,true, "in1");
		String expect = "GGTTACCTGCCGTGAGTAAA"; 
		//String iprimer = PCR.getUnmodifiedForwardPrimer(insert, genome);
		//System.out.println("Forward from insert oligo: " + iprimer);
		//System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.equals(iexpect)));
                PCR pcr = new PCR(genome);
                Primer p = pcr.getUnmodifiedForwardPrimer(insert);
                System.out.println("expected: " + expect);
                System.out.println("UnmodifiedFP insertion correct: " + expect.equals(p.seq));
                System.out.println("actual  : " + p.seq);
                
                Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
                p = pcr.getUnmodifiedForwardPrimer(mismatch);
                System.out.println("UnmodifiedFP mismatch correct: " + expect.equals(p.seq));
                System.out.println("actual  : " + p.seq);
                
                Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "del1");
                p = pcr.getUnmodifiedForwardPrimer(delete);
                System.out.println("UnmodifiedFP delete correct: " + expect.equals(p.seq));
                System.out.println("actual  : " + p.seq);
                

	}
	
	public static void testGetModifiedForwardPrimerFromOligo() throws Exception{
            String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
            PCR pcr = new PCR(genome);
            Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2,true, "in1");
            String iexpect = "GTTACCTGCCCCCCCAAAAA"; //base 90-99, then the insert
            Primer iprimer = pcr.getModifiedForwardPrimer(insert);
            System.out.println("Forward from insert oligo: " + iprimer.seq);
            System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.seq.equals(iexpect)));

            Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
            String mexpect = "GTTACCTGCCCCCCCAAAAA"; //80-89, then the new seq
            Primer mprimer = pcr.getModifiedForwardPrimer(mismatch);
            System.out.println("Forward from mismatch oligo: " + mprimer.seq);
            System.out.println("Test getForwardPrimerFromOligo mismatch : " + (mprimer.seq.equals(mexpect)));

            genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTGGGGGGGGGGAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
            pcr = new PCR(genome);
            Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "del1");
            String dexpect = "TTTTTTTTTTTTTTTAAAAA"; //bases 80-89,100-109
            Primer dprimer = pcr.getModifiedForwardPrimer(delete);
            System.out.println("Forward from deletion oligo: " + dprimer.seq);
            System.out.println("Test getForwardPrimerFromOligo deletion : " + (dprimer.seq.equals(dexpect)));

    }
	
	//replicate pcr primers used in http://www.sciencemag.org/content/suppl/2011/07/13/333.6040.348.DC1/Isaacs.SOM.pdf
	public static void testPCRbyReplicatingResults() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
                PCR pcr = new PCR(genome);
		Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730611, 4730612, 2, true, "TAG->TAA");
		String expectUnmodified = "AACTGGTTGTTAAGCAGTAG";
		String expectModified = "AACTGGTTGTTAAGCAGTAA";
		Primer unmodified = pcr.getUnmodifiedForwardPrimer(mismatch);
		Primer modified = pcr.getModifiedForwardPrimer(mismatch);
		System.out.println("Successfully replicated unmodified primer: " + unmodified.seq.equals(expectUnmodified));
		System.out.println("Successfully replicated modified primer: " + modified.seq.equals(expectModified));
	}
	
        public static void testGeneratePrimerSets() throws IOException, Exception{
            String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
            PCR pcr = new PCR(genome);
            ArrayList <Oligo> pool =  new ArrayList<Oligo>();
            pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "GACGATATGAAC", 7900, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "ATTTGACCT", 15426, 2,true, "oligo") );
            pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2,true, "oligo") );

            ArrayList<ArrayList<Primer>> res = pcr.generateAllPrimers(pool);
            System.out.print(res);
        }
        
        /*
	public static void testGetMASCPCRPrimersForOligo() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		Oligo insert = Oligo.InsertionFactory(genome, "TTTTTAAAAA", 100, 2,true, "in1");
		System.out.println("Generating primers. This may take a minute...");
		List<String> primers = PCR.getMASCPCRPrimers(insert,genome);
		System.out.println("Insertion Oligo: " + insert.getSequenceAsString());
		System.out.println("Insertion Primers:");
		System.out.println(primers);
		
		Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "del1");
		primers = PCR.getMASCPCRPrimers(delete,genome);
		System.out.println("Deletion Oligo: " + delete.getSequenceAsString());
		System.out.println("Deletion Primers:");
		System.out.println(primers);
		
	}
	*/

	
}
