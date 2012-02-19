package mage.optMage;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mage.Core.Oligo;
import utils.TextFile;

public class Parser {

	public static List<Integer> parse(String filepath, List<Oligo> pool) throws Exception {
		
		ArrayList<Integer> result = new ArrayList<Integer>(pool.size());
		List<String> optResults = parse(filepath);
		for (int ii = 0; ii< Math.min(pool.size(),optResults.size()) ;ii++) {
			result.add( adjust(optResults.get(ii), pool.get(ii) ) );
		}
		
		return result;
	}
	
	/**
	 * Identifies the oligo in the span that corresponds to the optMAGE selection
	 * for comparison
	 * 
	 * @param optmage
	 * @param ol
	 * @return
	 * @throws Exception
	 */
	public static int adjust(String optmage, Oligo ol) throws Exception{
		int shift = 0;
		int position = -1;
		
		// Find the oligo and return the position
		for (String poligo : ol.getPossibleOligos()){
			if (poligo.equalsIgnoreCase(optmage)) {position = shift; break;}
			shift ++; 
		}
		return position;
	}
	
	/**
	 * Helper function for parsing the optMAGE output
	 * 
	 * @param filepath		The File path of the optMAGE file
	 * @return				Returns a pool of oligos
	 * @throws IOException
	 */
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
