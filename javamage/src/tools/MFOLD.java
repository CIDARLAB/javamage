package tools;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class MFOLD {

	private List<String> list;
	private ArrayList<Double> results;

	public MFOLD (List<String> list) {
		setOligos(list);
		this.results = new ArrayList<Double>(list.size());
	}

	public List<String> getOligo() {
		return this.list;
	}

	public void setOligos(List<String> list) {
		this.list =  list;

	}

	public ArrayList<Double> run() {
				
		for ( String seq: list)
		{
			BufferedReader input;
			try{
				// Create and run mfold process
				ProcessBuilder pb = new ProcessBuilder(
						Constants.MFOLD,"--NA=DNA","--energyOnly","-q",seq.toString());
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

	public static void main(String[] args) {
		
		ArrayList<String> list = new ArrayList<String>();
		list.add(new String("atcgatcggctaggtaacagattaatctctcggagctgatacgac"));
		list.add(new String("TTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACGGGATTATTATTGGGCTCGAATCTACCGTCGATATTGCTGAGTCCACCC"));		
		MFOLD mf =  new MFOLD(list);
		
		for ( Double score : mf.run() ) {
			System.out.println(score);
		} 
	}

}
