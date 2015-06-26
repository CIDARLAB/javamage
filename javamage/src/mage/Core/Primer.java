package mage.Core;

/** Class to contain information and helper functions for a MASC PCR primer
 *
 * @author mquintin
 */
public class Primer {
    //is this the forward primer? If not, it's the reverse
    //the convention we're using is that the forward primers
    //are the ones that span the modified locus
    public boolean forward;
    //which direction is this primer set? if true, the pair follows the given
    //sequence orientation (5' -> 3'). If false, it's flipped (3' -> 5')
    public boolean sense;
    //sequence of this primer
    public String seq;
    //the Oligo this primer maps to
    public Oligo oligo;
    //the size of the amplicon created by this primer, if it's reverse
    public int amplicon;
    
    public Primer(String seq, Oligo oligo, int amplicon, boolean forward, boolean sense){
        this.forward = forward;
        this.sense = sense;
        this.seq = seq;
        this.oligo = oligo;
        this.amplicon = amplicon;
    }

    
}
