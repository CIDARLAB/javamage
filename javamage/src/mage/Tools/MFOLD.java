package mage.Tools;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


/**
 * A wrapper for mfold on a oligo and collecting the results 
 * @author Samir Ahmed
 *
 */
public class MFOLD {

	private List<String> list;
	private ArrayList<Double> results;

	/**
	 * Create a new mfold wrapper object
	 * @param list	A List of sequences to be run through mfold
	 */
	public MFOLD (List<String> list) {
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
	 * Sets the the list of sequences to run throuhg mfold
	 * 
	 * @param list An array list of sequences
	 */
	public void setOligos(List<String> list) {
		this.list =  list;

	}

	/**
	 * Will take the set of oligos and run them through Mfold
	 * 
	 * @returns a list of scores that correspond to the mfold results
	 */
	public ArrayList<Double> run() {
				
		for ( String seq: list)
		{
			BufferedReader input;
			try{
				// Create and run mfold process
				ProcessBuilder pb = new ProcessBuilder(
						Constants.MFOLD,"--NA=DNA","--energyOnly","-q",seq.toString());
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
		
		for ( Double score : mf.run() ) {
			//System.out.println(score);
		} 
	}

}
