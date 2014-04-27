package mage.Tools;

import java.util.ArrayList;
import java.util.List;

/**Class to support design of dsDNA primers
 * 
 * User will provide start and end positions, and the nucelotide sequence of a gene to be inserted at that position
 * 
 * Initial design will be to devise two 60 bp primers. The 40 5'-end bases will be homologous to the target region on the genome,
 * and the 20 3' bases will be homologous to the insertion region
 * 
 * See Daiguan Yu, Hilary Ellis et all, 2000, "An Efficient Recombination System for Chromosome Engineering 
 * in Escherichia coli", PNAS vol 97 no 11
 * 
 * See also Court lab recombineering FAQ: http://redrecombineering.ncifcrf.gov/Background.html
 * 
 * @author mquintin
 *
 */
public abstract class DSDNA {
	//TODO: this will break on deletions! allow for sequences less than 40bp
	static int genomeOverlap = 40;
	static int insertOverlap = 20;
	
	/**Return a List with two elements, for the two primers
	 * 
	 * @param genome genomic sequence
	 * @param sequence complete sequence to be inserted, though only the insertOverlap bases on either end actually matter
	 * @param leftpos start location of the replacement region on the genome. (first position in the genome ==1). This base will remain in the sequence
	 * @param rightpos location of the last base to be replaced
	 * @return
	 */
	public static List<String> getDSDNAPrimers(String genome, String sequence, int leftpos, int rightpos){
		//convert from 1-indexed to 0 indexed
		int left = leftpos -1;
		int right = rightpos -1;
		
		List<String> list = new ArrayList<String>();
		
		//if the insert sequence is blank or too short
		//currently displaying an error message- should look into how deletions/point mutations are generally done
		if (sequence.length() < 2 * insertOverlap){
			list.add("Error: The inserted sequence should be at least " + (2 * insertOverlap) + " bases long");
		}
		else{
			//insert sequence > 2*geneoverlap
			//build the left primer
			String genleft = genome.substring(left - genomeOverlap, left);
			String insleft = sequence.substring(0,insertOverlap);
			String lprimer = genleft + insleft;
			
			//build the right primer
			String genright = genome.substring(right,right + genomeOverlap);
			int len = sequence.length();
			String insright = sequence.substring(len - insertOverlap, len);
			String rprimer = SequenceTools.ReverseCompliment(genright + insright);
			
			list.add(lprimer);
			list.add(rprimer);
		}
		return list;
	}

}
