package test.Unit;

import java.util.List;

import test.Constants;
import mage.Core.Oligo;
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
		//testGetUnmodifiedForwardPrimerFromOligo();
		//testGetModifiedForwardPrimerFromOligo();
		//testGetMASCPCRPrimersForOligo();
		testPCRbyReplicatingResults();

	}
	
	public static void testGetUnmodifiedForwardPrimerFromOligo() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2,true, "in1");
		String iexpect = "TTCTGAACTGGTTACCTGCC"; 
		String iprimer = PCR.getUnmodifiedForwardPrimer(insert, genome);
		System.out.println("Forward from insert oligo: " + iprimer);
		System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.equals(iexpect)));

	}
	
	public static void testGetModifiedForwardPrimerFromOligo() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2,true, "in1");
		String iexpect = "GTTACCTGCCCCCCCAAAAA"; //base 90-99, then the insert
		String iprimer = PCR.getModifiedForwardPrimer(insert);
		System.out.println("Forward from insert oligo: " + iprimer);
		System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.equals(iexpect)));
		
		Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
		String mexpect = "GTTACCTGCCCCCCCAAAAA"; //80-89, then the new seq
		String mprimer = PCR.getModifiedForwardPrimer(mismatch);
		System.out.println("Forward from mismatch oligo: " + mprimer);
		System.out.println("Test getForwardPrimerFromOligo mismatch : " + (mprimer.equals(mexpect)));
		
		genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTGGGGGGGGGGAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
		Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "del1");
		String dexpect = "TTTTTTTTTTTTTTTAAAAA"; //bases 80-89,100-109
		String dprimer = PCR.getModifiedForwardPrimer(delete);
		System.out.println("Forward from deletion oligo: " + dprimer);
		System.out.println("Test getForwardPrimerFromOligo deletion : " + (dprimer.equals(dexpect)));
		
	}
	
	//replicate pcr primers used in http://www.sciencemag.org/content/suppl/2011/07/13/333.6040.348.DC1/Isaacs.SOM.pdf
	public static void testPCRbyReplicatingResults() throws Exception{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730611, 4730612, 2, true, "TAG->TAA");
		String expectUnmodified = "AACTGGTTGTTAAGCAGTAG";
		String expectModified = "AACTGGTTGTTAAGCAGTAA";
		String unmodified = PCR.getUnmodifiedForwardPrimer(mismatch, genome);
		String modified = PCR.getModifiedForwardPrimer(mismatch);
		System.out.println("Successfully replicated unmodified primer: " + unmodified.equals(expectUnmodified));
		System.out.println("Successfully replicated modified primer: " + modified.equals(expectModified));
	}
	
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
	

	
}
