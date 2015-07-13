package test;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import mage.Tools.OligoStats;

/**Compare ARE prediction methods to real-world data
 * See http://www.sciencemag.org/content/333/6040/348/suppl/DC1
 * 
 * @author Michael Quintin
 *
 */
public class TestReplacementPrediction {

	static ArrayList<Oligo> oligos;
	static String oligoFile = "C:/Users/mquintin/OneDrive/github/javamage/testing/validationOligos.txt";
	static String genomeDir = "C:/Users/mquintin/OneDrive/github/javamage/testing/";
	static String genomeFile = "MG1655.FASTA";
	
	static int oligoID = 1;
	
	public static void main(String[] args) throws Exception{
		populateOligos();
		//getStats();
		getSingleARE();
	}

	private static void populateOligos() throws Exception{
		
		
		/*BufferedReader reader = new BufferedReader( new FileReader(oligoFile));
		String line  = null;
		while( ( line = reader.readLine() ) != null ) {
			String[] arr = line.split("\t");
			String name = arr[0];
			String seq = arr[1];
			Oligo oligo = Oligo.MismatchFactory(genome, "TA", 1000, 1002, 1, true, name);
			
		}*/
		oligos = new ArrayList<Oligo>();
		for (int i = 0; i < 10; i++){
			oligos.add(getTestOligo());
		}
		
	}
	
	private static void getStats(){
		System.out.println(OligoStats.getDiscreteDiversityTable(oligos, 18));
	}
	
	/*all oligos in the validation set are 90bp, and contain a single base mismatch
	since the ARE function doesn't check sequence or position of the mismatch,
	we can repeat the same oligo
	*/
	
	private static Oligo getTestOligo() throws Exception{
		String genome = FASTA.readFFN(genomeDir, genomeFile);
		Oligo oligo = Oligo.MismatchFactory(genome, "T", 100, 101, 1, true, "oligo".concat(String.valueOf(oligoID)));
		oligoID++;
		return oligo;
	}
	
	private static void getSingleARE() throws Exception{
		Oligo oligo = getTestOligo();
		System.out.println(OligoStats.getARE(oligo));
	}
	
}
