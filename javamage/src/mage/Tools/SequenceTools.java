package mage.Tools;

import org.biojava3.core.sequence.DNASequence;
 
public class SequenceTools {
	
  public static String ReverseCompliment(String sequence) {
 
	  DNASequence seq = new DNASequence(sequence);
	  return seq.getReverseComplement().getSequenceAsString();
  }

}