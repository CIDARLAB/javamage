package mage.optMage;

import java.util.List;

import mage.Core.Oligo;
import utils.Plot;

public class Comparator {

	/**
	 * Plots the oligo for the given optMAGE positions
	 * @param ol				Oligo
	 * @param optMagePosition	optMAGE optimal position
	 */
	public static void plot(Oligo ol, int optMagePosition){
		
		// Set the oligos optmage position
		ol.setOptMagePosition(optMagePosition);
		
		// Create double arrays
		Double[] bgScores = ol.bgList().toArray(new Double[ol.bgList().size()]);
		Double[] dgScores = ol.dgList().toArray(new Double[ol.dgList().size()]);		
		Double[] boScores = ol.boList().toArray( new Double[ol.boList().size()]);
		
		// Create Plot object
		Plot pl = new Plot(1,2);
		
		// Get Line min and max points
		double line_min =(double) min(bgScores) - 1.0; 
		double line_max = (double) max(bgScores)+ 1.0;
		
		// Add BG plot and line for optMAGE position and MERLIN Position
		pl.addGraph( bgScores, "Genome Homology");
		pl.addLine(	new Double [] { (double) optMagePosition, (double) optMagePosition },  
					new Double [] { line_min , line_max }, 
					"optMAGE" );
		pl.addLine(	new Double [] { (double)ol.getOptimalPosition() , (double) ol.getOptimalPosition() }, 
					new Double [] { line_min , line_max },
					"MERLIN");
		
		// Add DG plot and lines for optMAGE position and MERLIN Position
		pl.addGraph(dgScores, "Free Energy");
		line_min =(double) min(dgScores) - .1; 
		line_max = (double) max(dgScores)+ .1;
		pl.addLine(	new Double [] { (double) optMagePosition, (double) optMagePosition },  
					new Double [] { line_min ,line_max  }, 
					"optMAGE" );
		pl.addLine(	new Double [] { (double)ol.getOptimalPosition() , (double) ol.getOptimalPosition() }, 
					new Double [] { line_min , line_max },
					"MERLIN");
		
		// Draw final plot
		pl.setToLines();
		pl.draw("Oligo_" + ol.getOligoId()+" DG_BG Comparison");

		pl = new Plot(1,1);
		
		// Add BO  plot and lines for optMAGE position and MERLIN Position
		pl.addGraph(boScores, "Oligo Pool Homology");
		line_min =(double) min(boScores) -1; 
		line_max = (double) max(boScores)+1;
		pl.addLine(	new Double [] { (double) optMagePosition, (double) optMagePosition },  
					new Double [] { line_min ,line_max  }, 
					"optMAGE" );
		pl.addLine(	new Double [] { (double)ol.getOptimalPosition() , (double) ol.getOptimalPosition() }, 
					new Double [] { line_min , line_max },
					"MERLIN");
		
		pl.setToLines();
		pl.draw("Oligo_" + ol.getOligoId()+" BO Comparison");
		
	}
	
	/**
	 * Returns the index of the minimum element of a given array
	 * @param darray 	Double Array
	 * @return			Integer index 
	 */
	private static double min(Double[] darray){
		double min = darray[0];
		for ( double dd : darray){
			min = Math.min(min,dd);
		}
		return min;
	}

	/**
	 * Returns the index of the mazimum element of a given array
	 * @param darray	Double Array
	 * @return			Integer Index
	 */
	private static double max(Double[] darray){
		double max = darray[0];
		for ( double dd : darray){
			max = Math.max(max,dd);
		}
		return max;
	}

	/**
	 * Given the filepath to the optmage OUToligo.txt file, this will load and plot the 
	 * difference between optMAGE and javaMAGE
	 * 
	 * @param filepath		the file path, name+directory
	 * @param pool			List of Oligos			
	 * @throws Exception	
	 */
	public static void compare(String filepath, List<Oligo> pool) throws Exception {  
		
		List<Integer> results = Parser.parse(filepath, pool);
		for ( int ii=0; ii< results.size() ; ii++) { 
			plot( pool.get(ii),results.get(ii));
		}
		
	}
	
	


}
