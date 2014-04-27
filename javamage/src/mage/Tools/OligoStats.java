package mage.Tools;

import mage.Core.Oligo;
import mage.Core.OligoType;

import java.lang.Math;
import java.util.ArrayList;


public class OligoStats{
	/*methods to build:
	*get ARE for one oligo
	*get aggregate ARE: sum of all AREs = # expected ARs per cycle
	*get aggregate ARE 2: probability of any ARs per cycle
	*iterate many cycles, find diversity
	*/
	
	/*Need:
	 * parameters that describe decay of ARE due to pool complexity
	 * 
	 * Aaron's document from last week with stats description
	 * 
	 */
	
	/**Get the Allelic Replacement Efficiency (ARE) for the given oligo.
	 * This is based on the similarity between the ssDNA oligo and the 
	 * corresponding primary sequence.
	 * 
	 * TODO: switch b to be the length of the change, not the span! oligos are assumed
	 * to be 90bp, which is already factored into the background probability
	 * 
	 * TODO: check if OPTMAGE formulae are the same thing we're using here
	 * optional: get formula for ARE based on free energy? [not as good a fit as this here]
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
	 * @param list List of Oligos
	 * @return
	 */
	public static double getAggregateAnyARE(ArrayList<Oligo> list){
		int n = list.size();
		double b = 1.59;
		double coef = Math.exp(b / n);
		System.out.println("Any coef: " + String.valueOf(coef));
		double p = 1;
		for (Oligo oligo : list){
			System.out.println("ARE for " + oligo.name + ": " + String.valueOf(getARE(oligo)));
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
		double subgamma = 1 - pn;
		double gamma = Math.pow(subgamma, (n-k));
		System.out.println("alpha,beta,gamma : " + alpha + ", "+beta +", " +gamma);
				
		//(n choose k) * (1-(1-ARE)^c)^k * (1-ARE)^[c(k-n)]
		//long alpha = comboChoose(n,k);
		//double bsub = Math.pow(1-are,c);
		//double beta = Math.pow(1-bsub, k); 
		//double beta = Math.pow(1 - Math.pow(1 - are,c),k);
		//double gamma = Math.pow(1 - are, (c * (k - n)));
		//System.out.println("n,k,ARE,c : "+n+", "+k+", "+are+", "+c);
		//System.out.println("alpha,beta,gamma : " + alpha + ", "+beta +", " +gamma);
		return (double)alpha * beta * gamma;
	}

	
	//same as above, done slightly different as a sanity check
	public static double lociProbability2(int k, int m, double are, int n){
		long alpha = comboChoose(k,m);
		double pn = 1 - Math.pow(1-are,n);
		double beta = Math.pow(pn,m);
		double subgamma = 1 - pn;
		double gamma = Math.pow(subgamma, (k-m));
		System.out.println("alpha2,beta2,gamma2 : " + alpha + ", "+beta +", " +gamma);
		return alpha * beta * gamma;
	}
	
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
	 public static String getDiversityTable(ArrayList<Oligo> oligos, int cycles){
		 String s = "";
		 double are = getAggregateAnyARE(oligos);
		 //double are = getAggregateSumARE(oligos);
		 int n = oligos.size();
		 for (int c = 1; c <= cycles; c++){
			 String row = "";
			 for (int i = 0; i <= n; i++){//i = k = the number of mutated loci
				 System.out.println("n,k,ARE,c: " + n + ","+i+","+are+","+c+"; "+lociProbability(n,i,are,c));
				 row = row.concat(" " + lociProbability(n,i,are,c));
			 }
			 s = s.concat(row + "\n");
		 }
		 return s;
	 }
	 
	//TODO:
	//fraction of loci modified in top clone
	//mean fraction of loci modified across whole population
	//output table results
	
	 /** Generate the graph image for the diversity trend
	  * 
	  * @param oligos
	  * @param cycles
	  */
	 public static void generateDiversityGraph(ArrayList<Oligo> oligos, int cycles){
		 String data = getDiversityTable(oligos, cycles);
		 generateDiversityGraph(data);
	 }
	 
	 /** Generate the graph image for the diversity trend
	  * 
	  * @param data
	  */
	 public static void generateDiversityGraph(String data){
		 
	 }
}