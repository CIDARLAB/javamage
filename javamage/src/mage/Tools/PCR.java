package mage.Tools;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mage.Core.Oligo;
import mage.Core.OligoType;
import mage.Core.Primer;

/**Class to generate MASC PCR primers for oligos
 * see http://wanglab.c2b2.columbia.edu/publications/2011_MIE_Wang.pdf
 * 
 * Assumptions
 * 	Begin the primer placing the first modified base at the 3' position
 * 	Direction doesn't matter- if the opposite strand gets better score, use that
 * 	Genome input is already 5' to 3'
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
 * Melting temperatures are determined by the melt.pl script.	
 * 
 * 
 * Confirm method for forward unmodified deletion primer: is currently only including 5 bases into the cut
 * This should be redone to put the leftmost deleted base at the end
 * 
 * @author mquintin
 *
 */

//TODO: Account for cases where the primers may extend past the start/end of the genome 
public abstract class PCR {
	private static List<Integer> ampliconLengths = Arrays.asList(100,150,200,250,300,400,500,700,850);
	private static int primerlength = 20; //generally 18-30
	//number of valid primers to explore, keeping the first modified base in the 3' half
	private static double maxshift = Math.floor(primerlength/2); 
	
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
	public static Primer getUnmodifiedForwardPrimer(Oligo oligo, String genome){
		//position of the last targeted base on the reference genome
		int end = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;
		
		//for deletions, try to cover as much of the deleted region as possible
		if (oligo.getOligoType() == OligoType.DELETION){
			end = oligo.getGenomeEnd() + oligo.target_position - (oligo.getLength() + oligo.getGenomeStart());
			end = oligo.target_position + oligo.getGenomeEnd() - oligo.getLength();
		}
		end = end - 1; //necessary to align this with the modified primers 
		String seq = genome.substring(end - primerlength, end);
                Primer primer = new Primer(seq, oligo, 0, true, true);
		return primer;
	}
        
        public static Primer getUnmodifiedAntisenseForwardPrimer(Oligo oligo, String genome){
		//position of the last targeted base on the reference genome
		int start = oligo.getGenomeStart() + oligo.target_position;
		
		String seq = genome.substring(start, primerlength + start);
                seq = SequenceTools.ReverseCompliment(seq);
                Primer primer = new Primer(seq, oligo, 0, true, false);
		return primer;
	}
	
	/**Get the forward PCR primer for the modified sequence, given the Oligo
	 * 
	 * @param genome
	 * @param start
	 * @return
	 */
	public static Primer getModifiedForwardPrimer(Oligo oligo){
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
		
		String seq = sequence.substring(start, start + primerlength);
                Primer primer = new Primer(seq, oligo, 0, true, true);
		return primer;
	}
        
        public static Primer getModifiedAntisenseForwardPrimer(Oligo oligo){
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
		
		String seq = sequence.substring(start, start + primerlength);
                Primer primer = new Primer(seq, oligo, 0, true, true);
		return primer;
	}

	
	/**Get a single PCR reverse primer
	 * 
	 * @param genome reference genome
	 * @param start start position for the forward primer on the reference genome
	 * @param amplength length of the amplicon, NOT the primer
	 * @return primer sequence as a String of length = primerlength
	 */
	/*public static Primer getDownstreamReversePrimer(String genome, Oligo oligo, int start, int amplength){
		String gen = genome;
		//gen = new StringBuffer(gen).reverse().toString();
		//int s = gen.length()-start; //where to start this primer
		int s = start + amplength;
		//if the primer is too close to the end, it will run off
		//take enough bases from the start, paste them onto the end
		if (gen.length() - s < amplength){
			gen = gen.concat(gen.substring(0,amplength));
		}
		String seq = gen.substring(s, s + primerlength);
                seq = SequenceTools.ReverseCompliment(seq);
                
                Primer primer = new Primer(seq, oligo, amplength, false, true);
                
		return primer;
	}
        
       	public static Primer getUpstreamReversePrimer(String genome, Oligo oligo, int start, int amplength){
		String gen = genome;
		//gen = new StringBuffer(gen).reverse().toString();
		//int s = gen.length()-start; //where to start this primer
		int s = start + oligo.target_length - amplength;
		//if the primer is too close to the end, it will run off
		//take enough bases from the start, paste them onto the end
		if (gen.length() - s < amplength){
			gen = gen.concat(gen.substring(0,amplength));
		}
		String seq = gen.substring(s, s + primerlength);
                seq = SequenceTools.ReverseCompliment(seq);
                
                Primer primer = new Primer(seq, oligo, amplength, false, true);
                
		return primer;
	}
*/

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
		//String ufp = getUnmodifiedForwardPrimer(oligo, genome);
		//primers.add(ufp);
		
		//String mfp = getModifiedForwardPrimer(oligo);
		//primers.add(mfp);

		//reverse
		//find the first modified base
		int start = oligo.getGenomeStart() + oligo.target_position;
		//deletions are a special case- put the gap 3/4 of the way into the sequence
		if (oligo.getOligoType().equals(OligoType.DELETION)){
			start = (int) (start - Math.floor(primerlength * 3 / 4));
		}

		for (int alen: ampliconLengths){ 
			//primers.add(getReversePrimer(genome,start, alen));
		}
		
		return primers;		
	}
        
        /**Generate a reverse primer
         * 
         * @param genome
         * @param oligo
         * @param len
         * @param baselength
         * @return 
         */
        public static Primer getSenseReversePrimer(String genome, Oligo oligo, int len, int baselength){
            //calculate the start location
            int start = oligo.getGenomeStart() + oligo.target_position;
            //deletions are a special case- put the gap 3/4 of the way into the sequence
            if (oligo.getOligoType().equals(OligoType.DELETION)){
                    start = (int) (start - Math.floor(primerlength * 3 / 4));
            }
            
            start += len; 
            String gen = genome;
            if (gen.length() < start + primerlength){
			gen = gen.concat(gen.substring(0,len + primerlength));
		}
            String seq = gen.substring(start, start + primerlength);
            seq = SequenceTools.ReverseCompliment(seq);

            Primer primer = new Primer(seq, oligo, baselength, false, true);

            return primer;
        }
        
        public static Primer getAntisenseReversePrimer(String genome, Oligo oligo, int len, int baselength){
            //calculate the start location
            int start = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;
            //deletions are a special case- put the gap 3/4 of the way into the sequence
            if (oligo.getOligoType().equals(OligoType.DELETION)){
                    start = (int) (start + Math.floor(primerlength * 3 / 4));
            }
            
            start -= len; 
            String gen = genome;
            if (0 > start - primerlength){
			gen = gen.concat(gen.substring(0,len + primerlength));
		}
            String seq = gen.substring(start, start + primerlength);
            
            Primer primer = new Primer(seq, oligo, baselength, false, false);

            return primer;
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
	public static List<List<String>> getMASCPCRPrimersV0(List<Oligo> pool) throws IOException{
		// First we read in the genome
		String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);
		
		List<List<String>> primerset = new ArrayList<List<String>>();
		for (Oligo oligo:pool){
			primerset.add(getMASCPCRPrimers(oligo,genome));
		}
		return primerset;
	}
	
        public static List<Primer> generateAllPrimers(List<Oligo> pool){
            List<Primer> list = new ArrayList();
            Primer primer = null;
            double range = 0.05; //how much variablility is acceptible in the position?
                //ex: if range = .05 the amplicon must be from 95% to 105% of the given length
            for (Oligo oligo : pool){
                //get the forward primers
                //5'->3' modified
                //5'->3' wildtype
                //3'->5' modified
                //3'->5' wildtype
                
                for (int amp : getAmpliconLengths()){
                    for (int i = (int) Math.ceil((1-range) * amp);
                            i <= Math.floor(amp * (1+range)); i++){
                        //get the 5'->3' set
                        
                        //get the 3'->5' set
                    }
                    
                }
            }
            return list;
        }


}
