package mage.Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mage.Core.Oligo;

/**Class to handle writing to files. This is the entry point for software (eg Magelet) using the 
 * Tools package to generate outputs via the DSDNA, OligoStats, and PCR classes
 * 
 * @author Michael Quintin
 *
 */
public class OutputTools {

	/**Generate the MASC PCR primers from the given oligos, and output them to the given file
	 * 
	 * @param pool
	 * @param dest
	 * @throws IOException 
	 */
	public static void generateMASCPCRFile(List<Oligo> pool, String dest) throws IOException{
		List<List<String>> primers = PCR.getMASCPCRPrimers(pool);
		List<String> names = new ArrayList<String>();
		for (Oligo oligo : pool){
			names.add(oligo.name);
		}
		
		writePrimerListsToFile(names, primers, dest);
	}
	
	/**Given the list of lists of MASC PCR primers, write to file
	 * 
	 * @param names list of names for each oligo
	 * @param primerLists list of lists of primers for each oligo
	 * @param dest file to create
	 * @throws IOException 
	 */
	public static void writePrimerListsToFile(List<String> names, List<List<String>> primerLists, String dest) throws IOException{
		if (names.size() != primerLists.size()){
			throw new IllegalArgumentException("Primer set names must agree with the number of primer sets");
		}
		
		File file = new File(dest);
		
		// if file doesn't exists, then create it
		if (!file.exists()) {
			file.createNewFile();
		}
		
		FileWriter fw = new FileWriter(file.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		//write header
		String header = ("Oligo\tUnmodifiedForward\tModifiedForward");
		for(int length : PCR.getAmpliconLengths()){
			header = header + "\t" + "Rev" + String.valueOf(length);
		}
		bw.write(header);
		
		//write each primer set
		for (int i = 0; i < primerLists.size(); i++){
			String str = names.get(i);
			for (String primer : primerLists.get(i)){
				str = str + "\t" + primer;
			}
			bw.newLine();
			bw.write(str);
		}
		bw.close();
	}
	
	/**Generate a file containing the two DSDNA primers for the given recombination. Note that the start and end
	 * are indexed with 1 (not 0) being the first base on the genome
	 * 
	 * An error message will be generated instead if the size of the replacement is too small (not at least twice the size of the
	 * overlap into the sequence to be inserted)
	 * 
	 * @param genome
	 * @param sequence Sequence of the insertion
	 * @param leftpos 1-indexed base position on the genome to begin the replacement
	 * @param rightpos 1-indexed base position on the genome to end the replacement
	 * @param dest
	 * @throws IOException 
	 */
	public static void generateDSDNAPrimerFile(String genome, String sequence, int leftpos, int rightpos, String dest) throws IOException{
		//make sure the start and end are in the right order
		int start = leftpos;
		int end = rightpos;
		if (leftpos > rightpos){
			start = rightpos;
			end = leftpos;
		}
		
		List<String> primers = DSDNA.getDSDNAPrimers(genome, sequence, start, end);
		writeDSDNAPrimersToFile(primers, dest);
	}
	
	/**Write the DSDNA primers to the given location
	 * 
	 * @param primers
	 * @param dest
	 * @throws IOException 
	 */
	public static void writeDSDNAPrimersToFile(List<String> primers, String dest) throws IOException{
		//if there are two elements in the primers list, it worked. otherwise it's an error message

		File file = new File(dest);
		
		// if file doesn't exists, then create it
		if (!file.exists()) {
			file.createNewFile();
		}
		
		FileWriter fw = new FileWriter(file.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		for (String str : primers){
			bw.write(str);
			bw.newLine();
		}
		bw.close();
	}
	
	public static void generateDiversityTrendTableFile(List<Oligo> oligos, int cycles, String dest) throws IOException{
		String table = OligoStats.getDiversityTable(oligos, cycles);
		System.out.println(table);
		writeDiversityTableToFile(table, dest);		
	}
	
	public static void writeDiversityTableToFile(String table, String dest) throws IOException{
		File file = new File(dest);
		
		// if file doesn't exists, then create it
		if (!file.exists()) {
			file.createNewFile();
		}
		
		FileWriter fw = new FileWriter(file.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		//table is already formatted
		bw.write(table);
		bw.close();
		
		
	}
}
