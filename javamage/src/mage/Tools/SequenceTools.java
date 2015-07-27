package mage.Tools;

import java.util.HashMap;
import org.biojava3.core.sequence.DNASequence;

/**
 * Static helper for reverse complimenting a sequence
 * @author Samir Ahmed
 *
 */
public class SequenceTools {
	
  public static String ReverseCompliment(String sequence) {
	  //remove whitespace
	  sequence = sequence.replaceAll("\\s", "");
	  
	  DNASequence seq = new DNASequence(sequence);
	  return seq.getReverseComplement().getSequenceAsString();
  }
  
  private static HashMap<Character,Character[]> baseCodes;
  
    /**
     * populate the map if it's not populated, and return it
     *
     * @return
     */
    public static HashMap<Character, Character[]> getIUPACCodes() {
        if (baseCodes == null) {
            baseCodes = new HashMap<>();
            baseCodes.put('A', new Character[]{'A', 'R', 'W', 'M', 'D', 'H', 'V', 'N'});
            baseCodes.put('C', new Character[]{'C', 'Y', 'S', 'M', 'B', 'H', 'V', 'N'});
            baseCodes.put('G', new Character[]{'G', 'R', 'S', 'K', 'B', 'D', 'V', 'N'});
            baseCodes.put('T', new Character[]{'T', 'Y', 'W', 'K', 'B', 'D', 'H', 'N'});
            baseCodes.put('U', new Character[]{'U', 'Y', 'W', 'K', 'B', 'D', 'H', 'N'});
            baseCodes.put('R', new Character[]{'R','D','V','N'});
            baseCodes.put('Y', new Character[]{'Y','B','H','N'});
            baseCodes.put('S', new Character[]{'S','B','V','N'});
            baseCodes.put('W', new Character[]{'W','D','H','N'});
            baseCodes.put('K', new Character[]{'K','B','D','N'});
            baseCodes.put('M', new Character[]{'M','H','V','N'});
            baseCodes.put('B', new Character[]{'B','N'});
            baseCodes.put('D', new Character[]{'D','N'});
            baseCodes.put('H', new Character[]{'H','N'});
            baseCodes.put('V', new Character[]{'V','N'});
            baseCodes.put('N', new Character[]{'N'});
            baseCodes.put('.', new Character[]{'.','-'});
            baseCodes.put('-', new Character[]{'.','-'});
        }
        return baseCodes;
    }

  /**Do the given chars representing a DNA base match according to IUPAC codes,
   * or is c1 in the set represented by c2?
   * Order dependent in terms of which base is the wildcard: baseMatch(A,R) 
   * returns true while (R,A) is false
   * @param c1
   * @param c2
   * @return 
   */
  public static boolean baseMatch(char c1, char c2){
      HashMap<Character,Character[]> map = getIUPACCodes();
      Character[] vals = map.get(Character.valueOf(c1));
      Character target = Character.valueOf(c2);
      boolean res = false;
      for (Character c : vals){
          if (target.equals(c)){
              res = true;
              break;
          }
      }
      return res;
  }



}