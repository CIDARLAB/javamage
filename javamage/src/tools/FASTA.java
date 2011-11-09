package tools;

import java.io.IOException;

import mage.Oligo;

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
	
	public static void writeFFN(String directory, String filename, String sequence) throws IOException {
		
		// Split the string by length 80
		String fasta = ">0\n"+sequence; 
		String [] seqSplit = fasta.split("(?<=\\G.{80})");
		StringBuilder sb = new StringBuilder();
		for (String ss :seqSplit) {	sb.append(ss); }
		
		TextFile.write(directory, filename, sb.toString());
	}
	
	public static void main(String[] args) throws IOException {
		
		System.out.println(readFFN(Constants.blastdirectory,Oligo.Genome).length() );
	}

}
