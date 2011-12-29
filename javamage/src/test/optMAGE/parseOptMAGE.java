package test.optMAGE;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import test.Constants;
import utils.TextFile;

public class parseOptMAGE {

	public static void main(String[] args) throws Exception {
	
		String genome = FASTA.readFFN(Constants.optMagedirectory, "genome.fasta");
		Oligo ol = Oligo.InsertionFactory(genome, "atcggg", 1000, 1, false);		
		
		List<String> shifts = parse(Constants.optMagedirectory+"test.txt");

		for (String ss : shifts){
					System.out.println("optMAGE Shift: " + ss + "\tAdjusted : "+ adjust(ss,ol));
		}
	}
	
	public static int adjust(String ol1, Oligo ol) throws Exception{
		int shift = 0;
		int position = -1;
		
		// Find the oligo and return the position
		for (String poligo : ol.getPossibleOligos()){
			if (poligo.equalsIgnoreCase(ol1)) {position = shift; break;}
			shift ++; 
		}
		return position;
	}
	
	public static List<String> parse(String filepath) throws IOException{
		
		String file = TextFile.read(filepath);
		
		String [] lines = file.split("\n");
		
		/*
		 * The output is in the form indicated below
		 * 
		 * ID	START	END		STRAND	REP	MUT	MSHIFT	dGss	OLIGOSIZE	MM_COUNT	INS_COUNT	DEL_COUNT	PRED_RE	OLIGO SEQ
		   pps	1785136	1785137	-		2	I	5		-9.4	90			0			1			0			11.83	ccgagttggttataccaaagcaccagcggtgacgagccattgttggacatxcgaacaatccttttgtgataaatgaacggtttgagaaac
		 */
		
		ArrayList<String> shifts = new ArrayList<String>();
		
		for (String line : lines) {
			if (!line.startsWith("ID") && !line.equals("")) {
				
				// Split by white space
				String [] values = line.split("\t");
				
				// Extract the shift amount
				shifts.add((values[values.length-1]));
			}
		}
		return shifts;
	}
}
