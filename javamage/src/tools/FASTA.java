package tools;

import java.io.IOException;

import utils.TextFile;

public class FASTA {

	
	/**
	 * Static function for reading, parsing and returning a string of nucleotides from a .ffn file
	 * 
	 * @param directory		Directory of FFN file
	 * @param filename		genome FFN file name
	 * @return				A String of nucleotides from the FFN FASTA File
	 * @throws IOException	
	 */
	public static String readFFN(String directory, String filename) throws IOException {
		
		// Load the String from file
		String rawFFN = TextFile.read(directory+filename);
		
		// Use string builder for increasing efficiency
		StringBuilder sb = new StringBuilder();
		String []  raw =  rawFFN.split("\n");
		for ( String ss : raw ) {
			if (!ss.startsWith(">") ){ sb.append(ss); }
		}
		
		// Generate new genome as String
		String genome = sb.toString();
		if (genome.isEmpty()) throw new IOException("[FASTA] Unable to parse genome");
		return genome;
	}
	
	public static void main(String[] args) throws IOException {
		
		System.out.println(readFFN(Constants.blastdirectory,"genome.ffn").length() );
	}

}
