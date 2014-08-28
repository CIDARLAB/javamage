package mage.Tools;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mage.Core.Oligo;
import mage.Core.OligoType;

/**Class to generate MASC PCR primers for oligos
 * see http://wanglab.c2b2.columbia.edu/publications/2011_MIE_Wang.pdf
 * 
 * Assumptions
 * 	Begin the primer placing the first modified base at the 3' position
 * 	Direction doesn't matter- if the opposite strand gets better score, use that
 * 	Genome input is already 3' to 5'
 * 	Uses a complete, circular genome
 * 	Primer for deletion: put the start 1/4 * primerlength before the deletion
 * 
 * Specific tasks to complete
 * 	Get primers in a specific direction at a specific length and strand for one oligo 
 * 	Get primers for one oligo (either strand/direction) at a specific length
 * 	Selection of which target has a primer of which length is arbitrary	
 * 
 * Special Considerations
 * 	amplicon lengths should be 100, 150, 200, 250, 300, 400, 500, 600, 700, 850 bps, 
 * 		which can produce clearly distinguishable bands on a 1.5% agarose gel
 * 		-How to handle more than 10 oligos?
 * 	
 * ignore free energy/melt temp FOR NOW	
 * 	What scoring can we use to call primers acceptable? (free energy, melt temp, etc)
 * 		How do you calculate melt temp? Free energy is obtained from MFOLD
 * 
 * 
 * Confirm method for forward unmodified deletion primer: is currently only including 5 bases into the cut
 * This should be redone to put the leftmost deleted base at the end
 * 
 * @author mquintin
 *
 */
public abstract class PCR {
	private static List<Integer> ampliconLengths = Arrays.asList(100,150,200,250,300,400,500,700,850);
	private static int primerlength = 20; //generally 18-30
		
	/**
	 * @return the ampliconLengths
	 */
	public static List<Integer> getAmpliconLengths() {
		return ampliconLengths;
	}

	
	/**Get the forward PCR primer for the UNMODIFIED sequence, 
	 * given a sequence and the start position
	 * 
	 * 
	 * @param oligo
	 * @param genome unmodified reference genome
	 * @return primer sequence as a String of length = primerlength
	 */
	public static String getUnmodifiedForwardPrimer(Oligo oligo, String genome){
		//position of the last targeted base on the reference genome
		int end = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;
		
		//for deletions, try to cover as much of the deleted region as possible
		if (oligo.getOligoType() == OligoType.DELETION){
			end = oligo.getGenomeEnd() + oligo.target_position - (oligo.getLength() + oligo.getGenomeStart());
			end = oligo.target_position + oligo.getGenomeEnd() - oligo.getLength();
		}
		end = end - 1; //necessary to align this with the modified primers 
		String primer = genome.substring(end - primerlength, end);
		return primer;
	}
	
	/**Get the forward PCR primer for the modified sequence, given the Oligo
	 * 
	 * @param genome
	 * @param start
	 * @return
	 */
	public static String getModifiedForwardPrimer(Oligo oligo){
		String sequence = oligo.getSequenceAsString();
		
		int start = primerlength; //just initializing this here, it will change 
						//3' end has the strongest effect on binding,
						//so that's what we track, then work backwards
						//default is just the start of the sequence
		
		//deletion is a special case- put the gap 3/4 of the way through the primer
		if (oligo.getOligoType().equals(OligoType.DELETION)){
			int gappos = oligo.target_position;
			start = (int) (gappos - Math.floor(primerlength * 3 / 4));
		}
		
		//find the last modified base
		else {
			start = oligo.target_position + oligo.target_length - primerlength;
			//this is the end of the target, which is fine for insertion. 
			//assuming the last indicated base is changed for mismatch
		}
		//primer start site
		//start = start + 1; //+1 b/c we'll want the base at index pos included
		
		String primer = sequence.substring(start, start + primerlength);
		return primer;
	}

	
	/**Get a single PCR reverse primer
	 * 
	 * @param genome reference genome
	 * @param start start position for the forward primer on the reference genome
	 * @param amplength length of the amplicon, NOT the primer
	 * @return primer sequence as a String of length = primerlength
	 */
	public static String getReversePrimer(String genome, int start, int amplength){
		String gen = genome;
		//gen = new StringBuffer(gen).reverse().toString();
		//int s = gen.length()-start; //where to start this primer
		int s = start + amplength;
		
		//System.out.println(gen);
		
		//if the primer is too close to the end, it will run off
		//take enough bases from the start, paste them onto the end
		if (gen.length() - s < amplength){
			gen = gen.concat(gen.substring(0,amplength));
		}
		
		String primer = gen.substring(s, s + primerlength);
		return SequenceTools.ReverseCompliment(primer);
	}

	/**Get the unmodified forward (first element), modified forward (second),
	 *  and reverse (third through eleventh elements) 
	 * PCR primers for this oligo
	 * 
	 * @param oligo
	 * @param genome
	 * @return
	 */
	public static List<String> getMASCPCRPrimers(Oligo oligo,String genome){
		List<String> primers = new ArrayList<String>();
		
		//forward
		String ufp = getUnmodifiedForwardPrimer(oligo, genome);
		primers.add(ufp);
		
		String mfp = getModifiedForwardPrimer(oligo);
		primers.add(mfp);

		//reverse
		//find the first modified base
		int start = oligo.getGenomeStart() + oligo.target_position;
		//deletions are a special case- put the gap 3/4 of the way into the sequence
		if (oligo.getOligoType().equals(OligoType.DELETION)){
			start = (int) (start - Math.floor(primerlength * 3 / 4));
		}

		for (int alen: ampliconLengths){ 
			primers.add(getReversePrimer(genome,start, alen));
		}
		
		return primers;		
	}

	/**Get all MASC PCR primers for the oligo pool as a list of lists. Each of the interior lists
	 * represents the primers for a single oligo, beginning with the forward unmodified and modified 
	 * primers and followed by reverse primers for each length in PCR.ampliconLengths
	 * 
	 * See Wang and Church 2011, Multiplexed Genome Engineering and Genotyping Methods: Applications
	 * for Synthetic Biology and Metabolic Engineering, Methods in Enzymology vol 498
	 * http://wanglab.c2b2.columbia.edu/publications/2011_MIE_Wang.pdf
	 * 
	 * @param pool List of Oligos
	 * @return PCR primers as List of Lists of Strings
	 * @throws IOException 
	 */
	public static List<List<String>> getMASCPCRPrimers(List<Oligo> pool) throws IOException{
		// First we read in the genome
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
		
		List<List<String>> primerset = new ArrayList<List<String>>();
		for (Oligo oligo:pool){
			primerset.add(getMASCPCRPrimers(oligo,genome));
		}
		return primerset;
	}
	


}
