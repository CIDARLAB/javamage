package mage.Tools;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


/**
 * A wrapper for melt.pl on a oligo and collecting the results 
 * See http://mfold.rna.albany.edu/?q=unafold-man-pages/melt.pl
 * @author Michael Quintin, adapted from code by Samir Ahmed
 *
 */
public class Melt {

	private List<String> list;
	private ArrayList<Double> results;

	//default parameters
	public static double tmin = 37; //temperature for energy minimization
	public static double na = 1; //molar sodium ion concentration
	public static double mg = 0; //molar magnesium ion concentration
	public static double ct = .0000002; //total molar strand concentration

	
	
	/**
	 * Create a new wrapper object
	 * @param list	A List of sequences to be run through the melt script
	 */
	public Melt (List<String> list) {
		setOligos(list);
		this.results = new ArrayList<Double>(list.size());
	}

	/**
	 * Returns the list of oligos
	 * @return
	 */
	public List<String> getOligo() {
		return this.list;
	}

	/**
	 * Sets the the list of sequences to run through melt
	 * 
	 * @param list An array list of sequences
	 */
	public void setOligos(List<String> list) {
		this.list =  list;

	}

	/**
	 * Will take the set of oligos and run them through the melt script
	 * 
	 * @returns a list of scores that correspond to the melting temperature results
	 */
	public ArrayList<Double> run() {
				
		for ( String seq: list)
		{
			BufferedReader input;
			try{
				// Create and run mfold process
				ProcessBuilder pb = new ProcessBuilder(
						//Constants.MFOLD,"--NA=DNA","--energyOnly","-q",seq.toString());
						Constants.melt,"--NA=DNA","--Ct="+String.valueOf(ct),"--sodium="+String.valueOf(na),
						"--magnesium="+String.valueOf(mg),seq.toString()); //TODO: inputting the sequence strings
				
				//System.out.println(pb.command());
				Process mfold = pb.start();

				// Extract result from stdout, and Store as a Double
				input = new BufferedReader(new InputStreamReader(mfold.getInputStream()));
				String dg = input.readLine();
				results.add(Double.parseDouble(dg));
				
				// Close the input stream
				input.close();
			}
			catch (Exception ee) {
				ee.printStackTrace();
				System.err.println("[MFOLD] Fatal Error Executing MFOLD");
			}
		}
		return this.results;
	}

	/** Test function
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		
		ArrayList<String> list = new ArrayList<String>();
		list.add(new String("atcgatcggctaggtaacagattaatctctcggagctgatacgac"));
		list.add(new String("TTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACGGGATTATTATTGGGCTCGAATCTACCGTCGATATTGCTGAGTCCACCC"));		
		MFOLD mf =  new MFOLD(list);
		
		/*for ( Double score : mf.run() ) {
			//System.out.println(score);
		} */
	}

	public static double getTmin() {
		return tmin;
	}

	public static void setTmin(double tmin) {
		Melt.tmin = tmin;
	}

	public static double getNa() {
		return na;
	}

	public static void setNa(double na) {
		Melt.na = na;
	}

	public static double getMg() {
		return mg;
	}

	public static void setMg(double mg) {
		Melt.mg = mg;
	}

	public static double getCt() {
		return ct;
	}

	public static void setCt(double ct) {
		Melt.ct = ct;
	}


}
