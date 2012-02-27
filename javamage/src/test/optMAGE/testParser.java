package test.optMAGE;

import java.util.List;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import mage.optMage.Parser;
import test.Constants;

public class testParser {

	public static void main(String[] args) throws Exception {
		
		String genome = FASTA.readFFN(Constants.optMagedirectory, "genome.fasta");
		Oligo ol = Oligo.InsertionFactory(genome, "atcggg", 1000, 1, false,"testoligo");		
		
		List<String> shifts = Parser.parse(Constants.optMagedirectory+"test.txt");

		for (String ss : shifts){
					System.out.println("optMAGE Shift: " + ss + "\tAdjusted : "+ Parser.adjust(ss,ol));
		}
	}
	
	
}
