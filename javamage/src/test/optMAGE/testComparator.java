package test.optMAGE;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import test.Constants;

public class testComparator {

	public static void main(String[] args) throws Exception {

		String genome = FASTA.readFFN(Constants.optMagedirectory, "genome.fasta");
		Oligo ol = Oligo.InsertionFactory(genome, "aggg", 1512, 1, false);	
		ol.calc_bg();
		ol.calc_dg();
		Comparator.plot(ol,27);
	}
	
}
