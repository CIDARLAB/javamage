package test.optMAGE;

import mage.Core.Oligo;
import mage.Tools.FASTA;
import test.Constants;
import utils.Plot;

public class compareBO_BG {


	public static void plot(Oligo ol, int optMagePosition){
		Double[] bgScores = ol.bgList().toArray(new Double[ol.bgList().size()]);
		Double[] dgScores  = ol.dgList().toArray(new Double[ol.dgList().size()]);		
		
		// Create Plot object
		Plot pl = new Plot(1,2);
		
		// Add plot and line
		pl.addGraph( bgScores, "Genome Homology");
		pl.addLine(	new Double [] { (double) optMagePosition, (double) optMagePosition },  
					new Double [] { (double) min(bgScores)-1 , (double) max(bgScores)+1  }, 
					"optMAGE" );
		
		// Add next plot and line
		pl.addGraph(dgScores, "Free Energy");
		pl.addLine(	new Double [] { (double) optMagePosition, (double) optMagePosition },  
					new Double [] { (double) min(dgScores)-1 , (double) max(dgScores)+1  }, 
					"optMAGE" );
		
		// Draw final plot
		pl.draw();
	}
	
	/**
	 * Returns the index of the minimum element of a given array
	 * @param darray 	Double Array
	 * @return			Integer index 
	 */
	private static int min(Double[] darray){
		double min=darray[0]; int count = 0; int index =0;
		for (Double dd: darray) {  
			if (min < dd) { 
				index =count; 
				min = dd; 
			} count++;
		}
		return index;
	}

	/**
	 * Returns the index of the mazimum element of a given array
	 * @param darray	Double Array
	 * @return			Integer Index
	 */
	private static int max(Double[] darray){
		double max=darray[0]; int count = 0; int index =0;
		for (Double dd: darray) { 
			if (max > dd) { 
				index =count; 
				max = dd; 
			} 
			count++; }
		return index;
	}

	
	public static void main(String[] args) throws Exception {

		String genome = FASTA.readFFN(Constants.optMagedirectory, "genome.fasta");
		Oligo ol = Oligo.InsertionFactory(genome, "atcggg", 1000, 1, false);	
		ol.calc_bg();
		ol.calc_dg();
		plot(ol,27);
	}

}
