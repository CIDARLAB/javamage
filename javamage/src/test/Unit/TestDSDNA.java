package test.Unit;

import java.io.IOException;
import java.util.List;

import test.Constants;
import mage.Tools.DSDNA;
import mage.Tools.FASTA;
import mage.Tools.SequenceTools;

/**A class to test dual-stranded DNA oligo generation
 * 
 * @author mquintin
 *
 */
public class TestDSDNA {
	
	static String gfp = "atgagcaaaggcgaagaactgtttaccggcgtggtgccgattctggtggaactggatggcgatgtgaacggccataaatttagcgtgagcggcgaaggcgaaggcgatgcgacctatggcaaactgaccctgaaatttatttgcaccaccggcaaactgccggtgccgtggccgaccctggtgaccacctttagctatggcgtgcagtgctttagccgctatccggatcatatgaaacagcatgatttttttaaaagcgcgatgccggaaggctatgtgcaggaacgcaccattttttttaaagatgatggcaactataaaacccgcgcggaagtgaaatttgaaggcgataccctggtgaaccgcattgaactgaaaggcattgattttaaagaagatggcaacattctgggccataaactggaatataactataacagccataacgtgtatattatggcggataaacagaaaaacggcattaaagtgaactttaaaattcgccataacattgaagatggcagcgtgcagctggcggatcattatcagcagaacaccccgattggcgatggcccggtgctgctgccggataaccattatctgagcacccagagcgcgctgagcaaagatccgaacgaaaaacgcgatcatatggtgctgctggaatttgtgaccgcggcgggcattacccatggcatggatgaactgtataaa";
	
	public static void main(String[] args) throws IOException{
		//testInsert();
		//testDelete();
		testTooShortError();
	}

	//insert GFP into a specific arbitrarily chosen location
	public static void testInsert() throws IOException{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		List<String> primers = DSDNA.getDSDNAPrimers(genome, gfp, 1000, 2000);
		String expectLeft = "GGACGCAACGGTTCCGACTACTCCGCGGCGGTGCTGGCTG" + "ATGAGCAAAGGCGAAGAACT";
		String expectRight = "GGTCATTGTTGACTGCACCTCCAGCCAGGCAGTGGCGGAT" + "GCATGGATGAACTGTATAAA";
		expectRight = SequenceTools.ReverseCompliment(expectRight);
		boolean leftCorrect = primers.get(0).toUpperCase().equals(expectLeft);
		boolean rightCorrect = primers.get(1).toUpperCase().equals(expectRight);
		System.out.println(primers.get(0));
		System.out.println("Left insertion primer correct : " + leftCorrect);
		System.out.println(primers.get(1));
		System.out.println("Right insertion primer correct : " + rightCorrect);
	}
	
	//delete 1000 bases- currently not implemented
	public static void testDelete() throws IOException{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		List<String> primers = DSDNA.getDSDNAPrimers(genome, "", 1000, 2000);
		String expectLeft = "GGACGCAACGGTTCCGACTACTCCGCGGCGGTGCTGGCTG";
		String expectRight = "GGTCATTGTTGACTGCACCTCCAGCCAGGCAGTGGCGGAT";
		expectRight = SequenceTools.ReverseCompliment(expectRight);
		boolean leftCorrect = primers.get(0).toUpperCase().equals(expectLeft);
		boolean rightCorrect = primers.get(1).toUpperCase().equals(expectRight);
		System.out.println(primers.get(0));
		System.out.println("Left deletion primer correct : " + leftCorrect);
		System.out.println(primers.get(1));
		System.out.println("Right deletion primer correct : " + rightCorrect);
		
	}
	
	//get an error message if the insert is not at least 40 bases
	public static void testTooShortError() throws IOException{
		String genome = FASTA.readFFN(Constants.blastdirectory,"ecoli.ffn");
		List<String> primers = DSDNA.getDSDNAPrimers(genome, "", 1000, 2000);
		System.out.println(primers);
	}
	
}
