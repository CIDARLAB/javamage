package test.optMAGE;

import utils.Plot;
import mage.Core.Oligo;

public class compareBO_BG {


	public static void plot(Oligo ol, int optMagePosition){
		Double[] bgScores = ol.bgList().toArray(new Double[ol.bgList().size()]);
		Double[] dgScores  = ol.dgList().toArray(new Double[ol.dgList().size()]);
		
		Plot pl = new Plot(1,2);
		
	}
	
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

	
	public static void main(String[] args) {

	}

}
