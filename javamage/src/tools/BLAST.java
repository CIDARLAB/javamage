package tools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import utils.TextFile;


/**
 * BLAST or Basic Local Search Alignment Tool, developed by Myers, Altschul, Gish, Miller and Lipman
 * is an algorithm for comparing primary biological sequence information, such as the amino-acid sequences 
 * of different proteins or the nucleotides of DNA sequences.
 * 
 * This Class is used to run Blast+ to blast nucleotides to determine their relative alignment score
 * 
 * @author Samir Ahmed
 *
 */
public class BLAST {

	private DataBase db;
	private String directory;
	private Boolean buildQuery;
	private String query;
	private String queryFileName = "query.fasta";
	private String subjectFFN;
	private ArrayList<Integer> index;
	/**
	 * 
	 * @param directory	String containing the directory in which the blast results are to be carried out
	 * @param subjectFFN	String containg the name of the .ffn file (i.e "genome.ffn") that will be used to make a sequence
	 */
	public BLAST(String directory, String subjectFFN){	

		this.subjectFFN = subjectFFN;
		this.query = "";
		this.buildQuery= false;
		this.db = new DataBase(subjectFFN);
		this.directory = directory;
		this.index = new ArrayList<Integer>();

		try{
			ProcessBuilder dbBuilder= new ProcessBuilder(makeBlastDBCommand());
			Process makeDB = dbBuilder.start(); 
			makeDB.waitFor();
		}
		catch (Exception ee){ee.printStackTrace();}

	}

	/**
	 * Given a collection of sequences, a set of queries are build for blasting against the subject
	 * @param queryList A linked list with each sequence to be queried
	 */
	public void setQuery(HashMap<Integer,String> queryMap)
	{	
		this.buildQuery = true;
		String queryTerms = "";
		//int count = 0;
		for( Entry<Integer, String> entry : queryMap.entrySet() )
		{
			queryTerms+="> Oligo_"+entry.getKey().toString();
			queryTerms+="\n";
			index.add(entry.getKey());
			// Split the string by length 80
			String [] seqSplit = entry.getValue().split("(?<=\\G.{80})");

			// Add each split string to the query term as a new line
			for (String singleLine : seqSplit){ queryTerms += singleLine + "\n";}
		}

		// Write this to file
		this.query = queryTerms;
		try {
			TextFile.write(this.directory,this.queryFileName,this.query);
		} catch (IOException e) { e.printStackTrace(); }
	}

	//	/**
	//	 * Set the query file name to a premade query file name
	//	 * 
	//	 * @param newQueryFileName .fasta file with the query sequence e.g "query.fasta" or "oligo1.fasta"
	//	 */
	//	public void setQuery(String newQueryFileName, ) {this.queryFileName = newQueryFileName;}

	/**
	 * Execute a Basic Local Alignment Search on the subject sequence with the query sequences 
	 * (made with the setQuery method)
	 * After calling the external process, this command will waitfor it to complete, parse and store the results
	 */
	public void run(){
		if (this.buildQuery){
			try{

				ProcessBuilder blastn = new ProcessBuilder(makeBlastNCommand());
				Process blast_standalone = blastn.start();
				blast_standalone.waitFor();

			}
			catch (Exception ee) {
				System.err.println("[BLAST] Fatal Error: Could Successfully Run blastn command"); 
				ee.printStackTrace(); 
			}

			try{
				String output = TextFile.read(this.directory+this.db.textFile);
				for( BlastResult br : parse(output) ) { System.out.println(br);
					
				}

			}
			catch (Exception ee) {
				System.err.println("[BLAST] Error Parsing BLAST Result");
				ee.printStackTrace();

			}
		}
	}

	/**
	 * Getter for the subject sequence that is being BLASTed against
	 * @return String with the subject sequence
	 */
	public String getSubject(){
		return this.db.name;
	}

	// /*
	public static void main(String[] args) {

		String directory =  "/Users/mockingbird/dropbox/research/optimization/blast/";
		String subjectFFN = "genome.ffn";

		HashMap <Integer,String> seqMap = new HashMap<Integer,String>();
		seqMap.put(12,"CAAATGTGCAGCATACGTATTTGCTCGGCGTGCTTGGTCTCTCGTACTTCTCCTGGAGATCAAGGAAATGTTTCTTGTCCAAGCGGACAG");
		seqMap.put(144,"AGCTCTCGATAGACGACTACGGGAAAATACGATCAGGACTCGGGACTACGATATACGACAGAAAGACGTACGATTTACGACTGCGCGCGA");

		BLAST bg = new BLAST(directory,subjectFFN);
		bg.setQuery(seqMap);
		bg.run();

	} 
	// */

	private List<BlastResult> parse(String output){
		String [] lines = output.split("\\n");
		ArrayList<BlastResult> results = new ArrayList<BlastResult>();
		Integer count = 0;

		for(String line: lines){

			if (line.startsWith("#")){
				if (line.contains("# Query:")) {
					count++;
					//System.out.println(count);
				}
			}
			else{
				String[] parameters = line.split("\\t");
				if (parameters.length == 7) {

					BlastResult br = new BlastResult( this.index.get(count-1) );
					br.qStart = Integer.parseInt(parameters[0]);
					br.qEnd = Integer.parseInt(parameters[1]);
					br.sStart = Integer.parseInt(parameters[2]);
					br.sEnd = Integer.parseInt(parameters[3]);
					br.match_sequence = parameters[4];
					br.evalue = Double.parseDouble(parameters[5]);
					br.bitscore = Double.parseDouble(parameters[6]);

					results.add(br);
				}
			}
		}
		System.out.println("[BLAST] Total Number of Queries = "+count);
		
		return results;
	}
	

	/**
	 * Helper Function for building the makeblastdb command:
	 * makeblastdb -in ./genome.ffn -dbtype nucl  -out test.db
	 * 
	 * @return A List containing the commands to execute the makeblastdb function
	 */
	private List<String> makeBlastDBCommand(){
		List<String> list= new ArrayList<String>();

		list.add(Constants.makeblastdb);
		list.add("-in");
		list.add(this.directory+this.subjectFFN);
		list.add("-dbtype");
		list.add("nucl");
		list.add("-out");
		list.add(this.directory+this.db.fileName);

		return list;
	}

	/**
	 * Helper Function for building the following command : 
	 * blastn -db ./test.db  -query ./oligo.fasta -out ./testout.txt -word_size 11 -evalue 10 -outfmt "7 ssequid ssac qstart qend sstart send qseq evalue bitscore"
	 * 
	 * @return A List of command arguments for the processbuilder
	 */
	private List<String> makeBlastNCommand(){
		List<String> list= new ArrayList<String>();

		// Building the following process command
		list.add(Constants.blastn);
		list.add("-db");
		list.add(directory+this.db.fileName);
		list.add("-query");
		list.add(directory+this.queryFileName);
		list.add("-out");
		list.add(directory+this.db.textFile);
		list.add("-word_size");
		list.add("11");
		list.add("-evalue");
		list.add("10");
		list.add("-outfmt");
		list.add("7 ssequid ssac qstart qend sstart send qseq evalue bitscore");

		return list;
	}

	/**
	 * 
	 * Simple Class for holding db information
	 * 
	 * @author Samir Ahmed
	 *
	 */
	private class DataBase {

		public String name;
		public String fileName;
		public String textFile;

		public DataBase(String subject) {
			this.name = subject.replaceAll("\\..*",""); 		// Remove file extension
			this.name = this.name.replaceAll(".fasta","");
			this.fileName = this.name + ".db";
			this.textFile = this.name + ".txt";
		}
	}

	private class BlastResult{

		public Integer oligoID;
		public String match_sequence;
		public Integer qStart;
		public Integer qEnd;
		public Integer sStart;
		public Integer sEnd;
		public Double  evalue;
		public Double bitscore;

		public BlastResult(Integer id) { 
			this.oligoID = id;
			this.qStart = -1;
			this.qEnd= -1;
			this.sStart = -1;
			this.sEnd = -1;
			this.match_sequence = "";
			this.bitscore = -1.0;
			this.evalue = -1.0;
		}
		
		@Override
		public String toString(){
			return oligoID.toString()+"\t"+qStart.toString()+"\t"+
					qEnd.toString()+"\t"+sStart.toString()+"\t"+
					sEnd.toString()+"\t"+match_sequence+"\t\t"+
					bitscore.toString()+"\t"+evalue.toString();
		}
		
	}
}
