package mage.Tools.Pcr;

import java.util.ArrayList;
import java.util.List;
import mage.Core.Oligo;
import mage.Core.OligoType;
import mage.Core.Primer;
import mage.Tools.SequenceTools;

/**
 * Class to handle PCR Primer creation and modification methods
 *
 * @author Mike Quintin
 */
public class PrimerFactory {

    static int defaultPrimerLength = 20; //generally 18-30
    String genome;
    int primerLength;

    //--------Forward Primers
    //convenience method to make a forward primer
    private static Primer forwardPrimer(String seq, Oligo oligo, int genome_start, boolean sense, boolean modified) {
        Primer p = new Primer(seq, oligo, 0, genome_start, true, sense, modified);
        return p;
    }

    /**
     * Get the forward PCR primer for the UNMODIFIED sequence, given a sequence
     * and the genome_start position
     *
     *
     * @param oligo
     * @param genome unmodified reference genome
     * @return primer sequence as a String of length = primerLength
     */
    public Primer getUnmodifiedForwardPrimer(Oligo oligo) {
        int genome_end = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;
        if (oligo.getOligoType() == OligoType.DELETION) {
            genome_end = oligo.getGenomeEnd() + oligo.target_position - (oligo.getLength() + oligo.getGenomeStart());
            genome_end = oligo.target_position + oligo.getGenomeEnd() - oligo.getLength();
        }
        genome_end = genome_end - 1;
        String seq = genome.substring(genome_end - primerLength, genome_end);
        Primer primer = forwardPrimer(seq, oligo, genome_end - primerLength, true, false);
        return primer;
    }

    /**
     * Get the forward PCR primer for the modified sequence, given the Oligo
     *
     * @param genome
     * @return
     */
    public Primer getModifiedForwardPrimer(Oligo oligo) {
        String sequence = oligo.getSequenceAsString();
        int relative_start = 0;
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            int gappos = oligo.target_position;
            relative_start = (int) (gappos - Math.floor(primerLength * 3 / 4));
        } else {
            relative_start = oligo.target_position + oligo.target_length - primerLength;
        }
        //handle the edge case of primer falling off the start of the oligo
        String preseq = "";
        int genomeStart = oligo.getGenomeStart() + relative_start;
        if (relative_start < 0){
            preseq = genome.substring(genomeStart, oligo.getGenomeStart());
            relative_start = 0;
        }
        sequence = preseq.concat(sequence);
        
        String seq = sequence.substring(relative_start, relative_start + primerLength);
        Primer primer = forwardPrimer(seq, oligo, genomeStart + relative_start, true, true);
        return primer;
    }
    
    /**
     * get the primer in the opposite direction. In this case, the last modified
     * base in the oligo will be the first base in the primer
     */
    public Primer getModifiedAntisenseForwardPrimer(Oligo oligo) {
        String sequence = oligo.getSequenceAsString();
        int relative_start = 0; //where the primer starts in relation to the oligo
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            int gappos = oligo.target_position;
            relative_start = (int) (gappos - Math.floor(primerLength * 1 / 4));
        } else {
            relative_start = oligo.target_position + oligo.target_length - primerLength;
        }
        //handle the edge case of primer falling off the start of the oligo
        String preseq = "";
        int genomeStart = oligo.getGenomeStart() + relative_start;
        if (relative_start < 0){
            preseq = genome.substring(genomeStart, oligo.getGenomeStart());
            relative_start = 0;
        }
        sequence = preseq.concat(sequence);
        
        String seq = sequence.substring(relative_start, relative_start + primerLength);
        seq = SequenceTools.ReverseCompliment(seq);
        Primer primer = forwardPrimer(seq, oligo, genomeStart, false, true);
        return primer;
    }

    public Primer getUnmodifiedAntisenseForwardPrimer(Oligo oligo) {
        int genome_start = oligo.getGenomeStart() + oligo.target_position;
        String seq = genome.substring(genome_start, primerLength + genome_start);
        seq = SequenceTools.ReverseCompliment(seq);
        Primer primer = forwardPrimer(seq, oligo, genome_start, false, false);
        return primer;
    }

    //---------Reverse Primers
    /**
     * Get a single PCR reverse primer on the sense strand
     *
     * @param oligo
     * @param genome_start start position for the forward primer on the
     * reference genome
     * @param amplength length of the amplicon, NOT the primer
     * @return primer sequence as a String of length = primerLength
     */
    public Primer getDownstreamReversePrimer(Oligo oligo, int genome_start, int amplength) {
        String gen = genome;
        int s = genome_start + amplength;
        if (gen.length() - s < amplength) {
            gen = gen.concat(gen.substring(0, amplength));
        }
        String seq = gen.substring(s, s + primerLength);
        seq = SequenceTools.ReverseCompliment(seq);
        Primer primer = new Primer(seq, oligo, amplength, s, false, true);
        return primer;
    }

    /**
     * deprecated: Generate a reverse primer
     *
     * @param genome
     * @param oligo
     * @param len
     * @param baselength
     * @return
     */
    /*public Primer getSenseReversePrimer(Oligo oligo, int len, int baselength) {
        int start = oligo.getGenomeStart() + oligo.target_position;
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            start = (int) (start - Math.floor(primerLength * 3 / 4));
        }
        start += len;
        String gen = genome;
        if (gen.length() < start + primerLength) {
            gen = gen.concat(gen.substring(0, len + primerLength));
        }
        String seq = gen.substring(start, start + primerLength);
        seq = SequenceTools.ReverseCompliment(seq);
        Primer primer = new Primer(seq, oligo, baselength, start, false, true);
        return primer;
    }*/

    public Primer getUpstreamReversePrimer(Oligo oligo, int genome_start, int amplength) {
        String gen = genome;
        int fp_relative_end = 0; //where the forward primer ends in relation to the oligo
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            int gappos = oligo.target_position;
            fp_relative_end = (int) (gappos + Math.floor(primerLength * 3 / 4));
        } else {
            fp_relative_end = oligo.target_position + oligo.target_length;
        }
        int s = genome_start + fp_relative_end - amplength;
        if (gen.length() - s < amplength) {
            gen = gen.concat(gen.substring(0, amplength));
        }
        if (s < 0) {
            gen = gen.substring(gen.length() + s) + gen.substring(0, amplength);
            s = 0;
        }
        String seq = gen.substring(s, s + primerLength);
        Primer primer = new Primer(seq, oligo, amplength, s, false, true);
        return primer;
    }

    //----------Utility Functions
    /**
     * extract the sequences from each primer and return their melting temps in
     * the same order
     *
     */
    public static Double[] getMeltingTemps(ArrayList<Primer> primers) {
        List<String> seqs = new ArrayList();
        for (Primer p : primers) {
            seqs.add(p.seq);
        }
        Double[] res = new Double[0];
        try {
            res = Melt.getMT(seqs);
        } catch (Exception ex) {
            System.err.println(ex.getMessage());
            System.err.println(ex.getStackTrace());
        }
        return res;
    }

    /**
     * Helper method to call the appropriate method for the Annealing tool which
     * requires being able to call for varying lengths and genome_start
     * positions
     *
     * @param oligo
     * @param amplicon
     * @param forward
     * @param sense
     * @param modified
     * @param shift how many bases to move the genome_start This trusts that
     * shift=0 for forward primers but doesn't check
     * @param length primer length
     * @return
     */
    public Primer getPrimer(Oligo oligo, int amplicon, boolean forward,
            boolean sense, boolean modified, int shift, int length) {
        primerLength = length;
        Primer primer = null;
        if (forward) {
            if (modified) {
                if (sense) {
                    primer = getModifiedForwardPrimer(oligo);
                } else {
                    primer = getModifiedAntisenseForwardPrimer(oligo);
                }
            } else {//unmodified
                if (sense) {
                    primer = getUnmodifiedForwardPrimer(oligo);
                } else {
                    primer = getUnmodifiedAntisenseForwardPrimer(oligo);
                }
            }
        } else { //reverse
            int start = oligo.target_position;
            int amplength = amplicon + shift;
            if (sense) {
                primer = getDownstreamReversePrimer(oligo, start, amplength);
            } else {
                start = start + oligo.target_length;
                primer = getUpstreamReversePrimer(oligo, start, amplength);
            }
            primer.amplicon = amplicon; //not the actual length- primers get sorted by this
        }
        return primer;
    }

    //-----------Constructors, Getters and Setters
    public PrimerFactory(String genome) {
        this.genome = genome;
        this.primerLength = defaultPrimerLength;
    }
    
    public PrimerFactory(String genome, int primerlength){
        this.genome = genome;
        this.primerLength = primerlength;
    }

    public static int getDefaultPrimerLength() {
        return defaultPrimerLength;
    }

    public static void setDefaultPrimerLength(int defaultPrimerLength) {
        PrimerFactory.defaultPrimerLength = defaultPrimerLength;
    }

    public int getPrimerLength() {
        return primerLength;
    }

    public void setGenome(String genome) {
        this.genome = genome;
    }

    public void setPrimerLength(int primerLength) {
        this.primerLength = primerLength;
    }

    public String getGenome() {
        return genome;
    }
}
