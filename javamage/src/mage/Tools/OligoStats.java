package mage.Tools;

import mage.Core.Oligo;
import mage.Core.OligoType;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
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
		double coef = Math.exp(-b / n);
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
	public static double getAggregateSumARE(List<Oligo> list){
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
	 
	 //########################################
	 //these methods are an alternate take on replacement prediction, and avoid using the 
	 //aggregate estimation
	 
	 /**find the probability that a cell carries a modification after the given number
	  * of cycles. Returns a list for the probability of each oligo in the input list, in
	  * the same order
	  * 
	  * @param pool
	  * @param cycle
	  * @return
	  */
	 public static ArrayList<Double> probabilitiesAfterCycles(List<Oligo> pool,int cycle){
		 ArrayList<Double> res = new ArrayList<Double>();
		 
		 int n = pool.size();
		 double b = 1.59; //empirically determined
		 double pf = Math.exp(-b / n);//pooling factor
		 
		 for (Oligo oligo:pool){
			 Double are = Double.valueOf(getARE(oligo));
			 are = are * pf; //include the pooling factor
			 Double val = 1-(Math.pow(1-are, cycle));
			 res.add(val);
		 }
		 return res;
	 }
	 
	 /**generate a table of probabilities that describes the odds of a cell having each oligo after
	  * each cycle up to the input maximum. Each outer list represents one cycle (from 1...cycle_max)
	  * each inner list is the probability for each oligo, in the same order they appear in the input list
	  *
 	  * In the original Python script for OptMAGE, this correlates to p_out
	  * @param pool
	  * @param cycle_max
	  * @return
	  */
	 public static List<List<Double>> probabilitiesForAllCycles(List<Oligo> pool, int cycle_max){
		 ArrayList<List<Double>> table = new ArrayList<List<Double>>();
		 for (int cycle = 1; cycle <= cycle_max; cycle++){
			 ArrayList<Double> row = probabilitiesAfterCycles(pool, cycle);
			 table.add(row);
		 }
		 return table;
	 }
	 
	 /**generate a table of probabilities for the number of transformations per cycle,
	  * or the Poisson binomial probability mass function
	  * These will be the height of discrete sections in a bar graph depicting the 
	  * diversification trend
	  * The outer list index is the number of completed cycles
	  * The inner list index is the number of transformations (from 0...nOligo)
	  * The value is the expected fraction of the population to have exactly that number
	  *  of transformations
	  * 
	  * 
	  * @param pool
	  * @param cycles
	  * @return
	  */
	 public static List<List<Double>> getDiscreteDiversityTable(List<Oligo> pool, int cycles){
		 List<List<Double>> p_table = probabilitiesForAllCycles(pool, cycles);
		 List<List<Double>> res = new ArrayList<List<Double>>();
		 for (List<Double> row : p_table){ ///row=p_list
			 //ArrayList<Double> resrow = new ArrayList<Double>();
			 Double[] resrow = new Double[pool.size()+1]; ///resrow=q_list
			 Arrays.fill(resrow, 0.0);
			 resrow[0] = 1.0;

			 int i = 0;
			 for (Double p : row){
				 int k = i + 1;
				 for (int j = 0; j < i + 1; j++){
					 resrow[k] = (1-p)*resrow[k] + p*resrow[k-1];
					 k -=1;
				 }
				 resrow[0] *= (1-p);
				 i += 1;
			 }
			 List<Double> list = Arrays.asList(resrow);
			 res.add(list);
		 }
		 return res;
	 }

	 /**generate a table of the cumulative distribution function for the number of
	  * transformations per cycle. These will be the step points in a bar graph showing
	  * the diversification trend 
	  * 
	  * @param pool
	  * @param cycles
	  * @return
	  */
	 public static List<List<Double>> getCumulativeDiversityTable(List<Oligo> pool, int cycles){
		 List<List<Double>> pdf_table = getDiscreteDiversityTable(pool, cycles);
		 List<List<Double>> res = new ArrayList<List<Double>>();
		 for (List<Double> row : pdf_table){
			 Double[] resrow = new Double[pool.size()+1];
			 Double val = 0.0;
			 for (int i = 0; i < resrow.length; i++){
				 val += row.get(i);
				 resrow[i] = val;
			 }
			 List<Double> list = Arrays.asList(resrow);
			 res.add(list);
		 }
		 return res;
	 }
	 
}