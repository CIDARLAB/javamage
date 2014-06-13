package mage.Tools;

import mage.Core.Oligo;
import mage.Core.OligoType;

import java.lang.Math;
import java.util.ArrayList;
import java.util.List;


public class OligoStats{

	/**Get the Allelic Replacement Efficiency (ARE) for the given oligo.
	 * This is based on the similarity between the ssDNA oligo and the 
	 * corresponding primary sequence.
	 * 
	 * @param oligo
	 * @return
	 */
	public static double getARE(Oligo oligo) throws RuntimeException{
		//See Wang and Church 2011, Optmage 0.9 Document
		//ARE = alpha * e^(beta(len-1))
		//where alpha and beta have been experimentally determined for
		//an oligo with a 90bp span and len is the length of the change
		
		OligoType type = oligo.getOligoType();
		double val = 0;
		int olen = oligo.target_length;

		//ARE formula depends on mutation type
		switch (type){
			case INSERTION:
				val = 0.15 * Math.exp(-0.075 * (olen-1));
				break;
			case DELETION:
				olen = Math.abs(oligo.getGenomeEnd() - oligo.getGenomeStart()) - oligo.span;
				val = 0.23 * Math.exp(-0.058 * (olen-1));
				break;
			case MISMATCH:
				val = 0.26 * Math.exp(-0.135 * (olen-1));
				break;
			default:
				throw new RuntimeException("Oligo does not have a legal OligoType");
		}
		return val;
	}
	
	/**Get the number of cycles necessary for the given oligo to be expressed at a given
	 * frequency in the population
	 * 
	 * @param oligo
	 * @param freq The expression frequency, between 0 and 1
	 * @return
	 */
	public static double getCyclesForFrequency(Oligo oligo, double freq){
		if (freq < 0 || freq > 1){
			throw new IllegalArgumentException("Replacement frequency must be between 0 and 1");
		}
		
		double are = getARE(oligo);
		double n = Math.log(1 - freq) / Math.log(1 - are);
		return  n;
	}
	
	/**Get the aggregate ARE, the probability of any mutations occurring this cycle,
	 * calculated as 1 minus the product of the probabilities for each oligo, 
	 * which are multiplied by the pooling factor to account for inter-oligo interference 
	 * 
	 * @param oligos List of Oligos
	 * @return
	 */
	public static double getAggregateAnyARE(List<Oligo> oligos){
		int n = oligos.size();
		double b = 1.59; //empirically determined
		double coef = Math.exp(b / n);
		//System.out.println("Any coef: " + String.valueOf(coef));
		double p = 1;
		for (Oligo oligo : oligos){
			//System.out.println("ARE for " + oligo.name + ": " + String.valueOf(getARE(oligo)));
			p = p * (1- (coef * getARE(oligo)));
		}
		return 1 - (p);
	}
	
	/**[Whats the simple explanation for what this means?]
	 * [Do both aggregate AREs use the pool depletion formula?]
	 * 
	 * @param list List of Oligos
	 * @return
	 */
	public static double getAggregateSumARE(ArrayList<Oligo> list){
		int n = list.size();
		double b = 1.59; //empirically determined
		double coef = Math.exp(b / n);
		double p = 0;
		for (Oligo oligo : list){
			p = p + (coef * getARE(oligo));
		}
		return p;
	}
	
	/**Probability that a clone has been modified at k of n loci
	 * after c cycles with the given aggregate ARE
	 * 
	 * P = (n choose k) * (1-(1-ARE)^c)^k * (1-ARE)^(c(n-k))
	 * 
	 * @param n total loci = number of oligos
	 * @param k  number of mutated loci
	 * @param are allelic replacement efficiency for this locus (or the pool)
	 * @param c cycle number
	 * @return
	 */
	public static double lociProbability(int n, int k, double are, int c){
		long alpha = comboChoose(n,k);
		double pn = 1 - Math.pow(1-are,c);
		double beta = Math.pow(pn,k);
		//double subgamma = 1 - are;
		double gamma = Math.pow((1 - are), (c * (n-k)));
		//System.out.println("n,k,are,c : " + n + ", " + k + ", " + are + ", " + c);
		//System.out.println("alpha,beta,gamma : " + alpha + ", "+beta +", " +gamma);
				
		return (double)alpha * beta * gamma;
	}

	
	//same as above, done slightly different as a sanity check
	/*public static double lociProbability2(int k, int m, double are, int n){
		long alpha = comboChoose(k,m);
		double pn = 1 - Math.pow(1-are,n);
		double beta = Math.pow(pn,m);
		double subgamma = 1 - pn;
		double gamma = Math.pow(subgamma, (k-m));
		System.out.println("alpha2,beta2,gamma2 : " + alpha + ", "+beta +", " +gamma);
		return alpha * beta * gamma;
	}*/
	
	//the next two functions implement n choose k
	//from http://professorjava.weebly.com/factorialchoose.html
	private static long factorial(int a){ //non-recursive factorial, returns long
	    int answer=1; 
	    for(int i=1;i<=a;i++){ //starts from 1, multiplies up, such as 1*2*3...
	      answer*=i;
	    }
	    return(answer); 
	  }
	 private static long comboChoose(int n, int k){ //combinatorics function. takes n,k
	    return(factorial(n)/(factorial(k)*factorial(n-k))); //definition of nCk.
	  }
	
	 /**Get a table of percent of the population mutated at 0-n loci after 1-c cycles,
	  * where n is the number of oligos in the input list
	  * 
	  * Uses the multiplicative concept of pooled ARE
	  * 
	  * @param oligos
	  * @param cycles
	  * @return
	  */
	 public static String getDiversityTable(List<Oligo> oligos, int cycles){
		 String s = "";
		 double are = getAggregateAnyARE(oligos);
		 //System.out.println("Aggregate ARE : " + are);
		 //double are = getAggregateSumARE(oligos);
		 int n = oligos.size();
		 for (int c = 1; c <= cycles; c++){
			 String row = "";
			 for (int i = 0; i <= n; i++){//i = k = the number of mutated loci
				 //System.out.println("n,k,ARE,c: " + n + ","+i+","+are+","+c+"; "+lociProbability(n,i,are,c));
				 row = row.concat(lociProbability(n,i,are,c) + "\t");
			 }
			 s = s.concat(row.trim() + System.getProperty("line.separator"));
		 }
		 return s;
	 }

	 //moved to OutputTools
	 /*/** Generate the graph image for the diversity trend
	  * 
	  * @param oligos
	  * @param cycles
	  */
	 /*public static void generateDiversityGraph(ArrayList<Oligo> oligos, int cycles){
		 String data = getDiversityTable(oligos, cycles);
		 generateDiversityGraph(data);
	 }
	 
	 /** Generate the graph image for the diversity trend
	  * 
	  * @param data
	  */
	 /*public static void generateDiversityGraph(String data){
		 
	 }*/
}