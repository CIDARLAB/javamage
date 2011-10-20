package tools;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.sequence.DNASequence;


public class MFOLD {

	private List<DNASequence> list;

	public MFOLD (List<DNASequence> list) {
		setOligos(list);
	}

	public List<DNASequence> getOligo() {
		return this.list;
	}

	public void setOligos(List<DNASequence> list) {
		this.list =  list;

	}

	public void run() {

		for ( DNASequence seq: list)
		{
			BufferedReader input;
			try{
				// Create and run mfold process
				ProcessBuilder pb = new ProcessBuilder(
						Constants.MFOLD,"--NA=DNA","--energyOnly","-q",seq.toString());
				Process mfold = pb.start();

				// Extract result from stdout
				input = new BufferedReader(new InputStreamReader(mfold.getInputStream()));
				String result = input.readLine();
				System.out.println(result);
				//Double.parseAsDouble(result));
				// Close the input stream
				input.close();
			}
			catch (Exception ee) {
				ee.printStackTrace();
				System.err.println("[MFOLD] Fatal Error Executing MFOLD");
			}
		}
	}

	public static void main(String[] args) {
		ArrayList<DNASequence> list = new ArrayList<DNASequence>();
		list.add(new DNASequence("atcgatcggctaggtaacagattaatctctcggagctgatacgac"));
		list.add(new DNASequence("ttggttataccaaagcaccagcggtgacgagccattgttggacatcgaacaatccttttgtgataaatgaacggtttgagaaacacatt"));		
		MFOLD mf =  new MFOLD(list);
		mf.run();
	}

}
