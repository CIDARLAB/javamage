package mage.Tools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.DoubleStream;
import java.util.stream.DoubleStream.Builder;

import mage.Core.Oligo;
import mage.Core.OligoType;
import mage.Core.Primer;

/**
 * Class to generate MASC PCR primers for oligos see
 * http://wanglab.c2b2.columbia.edu/publications/2011_MIE_Wang.pdf
 *
 * Assumptions Begin the primer placing the first modified base at the 3'
 * position Direction doesn't matter- if the opposite strand gets better score,
 * use that Genome input is already 5' to 3' Uses a complete, circular genome
 * Primer for deletion: put the start 1/4 * primerlength before the deletion
 *
 * Specific tasks to complete Get primers in a specific direction at a specific
 * length and strand for one oligo Get primers for one oligo (either
 * strand/direction) at a specific length Selection of which target has a primer
 * of which length is arbitrary
 *
 * Special Considerations amplicon lengths should be 100, 150, 200, 250, 300,
 * 400, 500, 600, 700, 850 bps, which can produce clearly distinguishable bands
 * on a 1.5% agarose gel -How to handle more than 10 oligos?
 *
 * Melting temperatures are determined by the melt.pl script.
 *
 *
 * TODO: Confirm method for forward unmodified deletion primer: is currently
 * only including 5 bases into the cut This should be redone to put the leftmost
 * deleted base at the end
 *
 * TODO: Confirm method for modified antisense forward primer. If the 3' buffer
 * length is shorter than the primer length, primer may run off the oligo span.
 * Need to add part of the genome seq to cover this case.
 *
 * @author mquintin
 *
 */
//TODO: Account for cases where the primers may extend past the start/end of the genome 
public class PCR {

    private String genome;
    private int primerlength;
    private static List<Integer> ampliconLengths = Arrays.asList(100, 150, 200, 250, 300, 400, 500, 700, 850);
    private static int defaultPrimerLength = 20; //generally 18-30
    //legal range of primer lengths
    private static double maxshift = Math.floor(defaultPrimerLength * 0.2);
    private static double mtrange = 1.5;//how far from the final temp is OK?

    public static boolean forceTemp = true; //if true, primers are based on the targetTemp
    //if false, primers are based on the average of the forward primer temps
    public static double targetTemp = 60.0;//default target melting temp in degrees C

    public static boolean getForceTemp() {
        return forceTemp;
    }

    public static void setForceTemp(boolean forceTemp) {
        PCR.forceTemp = forceTemp;
    }

    public static double getTargetTemp() {
        return targetTemp;
    }

    public static void setTargetTemp(double targetTemp) {
        PCR.targetTemp = targetTemp;
    }

    public PCR(String genome) {
        this.genome = genome;
        this.primerlength = defaultPrimerLength;
    }

    /**
     * @return the ampliconLengths
     */
    public static List<Integer> getAmpliconLengths() {
        return ampliconLengths;
    }

    /**
     * Get the forward PCR primer for the UNMODIFIED sequence, given a sequence
     * and the start position
     *
     *
     * @param oligo
     * @param genome unmodified reference genome
     * @return primer sequence as a String of length = primerlength
     */
    public Primer getUnmodifiedForwardPrimer(Oligo oligo) {
        //position of the last targeted base on the reference genome
        int end = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;

        //for deletions, try to cover as much of the deleted region as possible
        if (oligo.getOligoType() == OligoType.DELETION) {
            end = oligo.getGenomeEnd() + oligo.target_position - (oligo.getLength() + oligo.getGenomeStart());
            end = oligo.target_position + oligo.getGenomeEnd() - oligo.getLength();
        }
        end = end - 1; //necessary to align this with the modified primers 
        String seq = genome.substring(end - primerlength, end);
        Primer primer = forwardPrimer(seq, oligo, end - primerlength, true, false);
        return primer;
    }

    public Primer getUnmodifiedAntisenseForwardPrimer(Oligo oligo) {
        //position of the last targeted base on the reference genome
        int start = oligo.getGenomeStart() + oligo.target_position;

        String seq = genome.substring(start, primerlength + start);
        seq = SequenceTools.ReverseCompliment(seq);
        Primer primer = forwardPrimer(seq, oligo, start, false, false);
        return primer;
    }

    /**
     * Get the forward PCR primer for the modified sequence, given the Oligo
     *
     * @param genome
     * @param start
     * @return
     */
    public Primer getModifiedForwardPrimer(Oligo oligo) {
        String sequence = oligo.getSequenceAsString();

        int start = primerlength; //just initializing this here, it will change 
        //3' end has the strongest effect on binding,
        //so that's what we track, then work backwards
        //default is just the start of the sequence

        //deletion is a special case- put the gap 3/4 of the way through the primer
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            int gappos = oligo.target_position;
            start = (int) (gappos - Math.floor(primerlength * 3 / 4));
        } //find the last modified base
        else {
            start = oligo.target_position + oligo.target_length - primerlength;
                    //this is the end of the target, which is fine for insertion. 
            //assuming the last indicated base is changed for mismatch
        }
            //primer start site
        //start = start + 1; //+1 b/c we'll want the base at index pos included

        String seq = sequence.substring(start, start + primerlength);
        Primer primer = forwardPrimer(seq, oligo, start, true, true);
        return primer;
    }

    /**
     * get the primer in the opposite direction. In this case, the last modified
     * base in the oligo will be the first base in the primer
     */
    public Primer getModifiedAntisenseForwardPrimer(Oligo oligo) {
        String sequence = oligo.getSequenceAsString();

        int start = primerlength; //just initializing this here, it will change 
        //3' end has the strongest effect on binding,
        //so that's what we track, then work backwards
        //default is just the start of the sequence

        //deletion is a special case- put the gap 3/4 of the way through the primer
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            int gappos = oligo.target_position;
            start = (int) (gappos - Math.floor(primerlength * 3 / 4));
        } //find the last modified base
        else {
            start = oligo.target_position + oligo.target_length - primerlength;
                    //this is the end of the target, which is fine for insertion. 
            //assuming the last indicated base is changed for mismatch
        }
            //primer start site
        //start = start + 1; //+1 b/c we'll want the base at index pos included

        String seq = sequence.substring(start, start + primerlength);
        seq = SequenceTools.ReverseCompliment(seq);
        Primer primer = forwardPrimer(seq, oligo, start, false, true);
        return primer;
    }

    /**
     * Get a single PCR reverse primer
     *
     * @param oligo
     * @param start start position for the forward primer on the reference
     * genome
     * @param amplength length of the amplicon, NOT the primer
     * @return primer sequence as a String of length = primerlength
     */
    public Primer getDownstreamReversePrimer(Oligo oligo, int start, int amplength) {
        String gen = genome;
            //gen = new StringBuffer(gen).reverse().toString();
        //int s = gen.length()-start; //where to start this primer
        int s = start + amplength;
            //if the primer is too close to the end, it will run off
        //take enough bases from the start, paste them onto the end
        if (gen.length() - s < amplength) {
            gen = gen.concat(gen.substring(0, amplength));
        }
        String seq = gen.substring(s, s + primerlength);
        seq = SequenceTools.ReverseCompliment(seq);

        Primer primer = new Primer(seq, oligo, amplength, start, false, true);

        return primer;
    }

    public Primer getUpstreamReversePrimer(Oligo oligo, int start, int amplength) {
        String gen = genome;
            //gen = new StringBuffer(gen).reverse().toString();
        //int s = gen.length()-start; //where to start this primer
        int s = start + oligo.target_length - amplength;
            //if the primer is too close to the end, it will run off
        //take enough bases from the start, paste them onto the end
        if (gen.length() - s < amplength) {
            gen = gen.concat(gen.substring(0, amplength));
        }
        if (s < 0) { //handle spanning the origin
            gen = gen.substring(gen.length() + s) + gen.substring(0, amplength);
            s = 0;
        }
        String seq = gen.substring(s, s + primerlength);
        //seq = SequenceTools.ReverseCompliment(seq);

        Primer primer = new Primer(seq, oligo, amplength, start, false, true);

        return primer;
    }

    /**
     * Get the unmodified forward (first element), modified forward (second),
     * and reverse (third through eleventh elements) PCR primers for this oligo
     *
     * @param oligo
     * @param genome
     * @return
     */
    public List<String> getMASCPCRPrimers(Oligo oligo, String genome) {
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
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            start = (int) (start - Math.floor(primerlength * 3 / 4));
        }

        for (int alen : ampliconLengths) {
            //primers.add(getReversePrimer(genome,start, alen));
        }

        return primers;
    }

    /**
     * Generate a reverse primer
     *
     * @param genome
     * @param oligo
     * @param len
     * @param baselength
     * @return
     */
    public Primer getSenseReversePrimer(Oligo oligo, int len, int baselength) {
        //calculate the start location
        int start = oligo.getGenomeStart() + oligo.target_position;
        //deletions are a special case- put the gap 3/4 of the way into the sequence
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            start = (int) (start - Math.floor(primerlength * 3 / 4));
        }

        start += len;
        String gen = genome;
        if (gen.length() < start + primerlength) {
            gen = gen.concat(gen.substring(0, len + primerlength));
        }
        String seq = gen.substring(start, start + primerlength);
        seq = SequenceTools.ReverseCompliment(seq);

        Primer primer = new Primer(seq, oligo, baselength, start, false, true);

        return primer;
    }

    public Primer getAntisenseReversePrimer(Oligo oligo, int len, int baselength) {
        //calculate the start location
        int start = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;
        //deletions are a special case- put the gap 3/4 of the way into the sequence
        if (oligo.getOligoType().equals(OligoType.DELETION)) {
            start = (int) (start + Math.floor(primerlength * 3 / 4));
        }

        start -= len;
        String gen = genome;
        if (0 > start - primerlength) {
            gen = gen.concat(gen.substring(0, len + primerlength));
        }
        String seq = gen.substring(start, start + primerlength);

        Primer primer = new Primer(seq, oligo, baselength, start, false, false);

        return primer;
    }

    /**
     * Get all MASC PCR primers for the oligo pool as a list of lists. Each of
     * the interior lists represents the primers for a single oligo, beginning
     * with the forward unmodified and modified primers and followed by reverse
     * primers for each length in PCR.ampliconLengths
     *
     * See Wang and Church 2011, Multiplexed Genome Engineering and Genotyping
     * Methods: Applications for Synthetic Biology and Metabolic Engineering,
     * Methods in Enzymology vol 498
     * http://wanglab.c2b2.columbia.edu/publications/2011_MIE_Wang.pdf
     *
     * @param pool List of Oligos
     * @return PCR primers as List of Lists of Strings
     * @throws IOException
     */
    public List<List<String>> getMASCPCRPrimersV0(List<Oligo> pool) throws IOException {
        // First we read in the genome
        String genome = FASTA.readFFN(Oligo.Directory, Oligo.Genome);

        List<List<String>> primerset = new ArrayList<List<String>>();
        for (Oligo oligo : pool) {
            primerset.add(getMASCPCRPrimers(oligo, genome));
        }
        return primerset;
    }

    /**
     * get a list of lists of primers. The first list is the forward primers,
     * each subsequent list is a set of reverse primers that is grouped together
     * for one MASC PCR reaction. The forward primers should be paired with
     * these sets by checking the Primer.oligo field
     *
     * @param pool a List of Oligos
     * @return
     * @throws IOException
     */
    public ArrayList<ArrayList<Primer>> generateAllPrimers(List<Oligo> pool) throws IOException {
        ArrayList<Primer> reversePrimers = new ArrayList();

        ArrayList<Primer> acceptedForward = generateForwardSet(pool);

        //get the reverse primers that correspond with the selected forwards
        for (int amp : ampliconLengths) {
            for (Primer p : acceptedForward) {
                Oligo o = p.oligo;
                if (p.sense) {
                    reversePrimers.add(getDownstreamReversePrimer(o, p.start, amp));
                } else {
                    reversePrimers.add(getUpstreamReversePrimer(o, p.start, amp));
                }
            }
        }
        ArrayList<Primer> acceptedReverse = correctOutlyingPrimers(reversePrimers,true);

        //remove any reverse primer that is still unacceptable
        /*for(Primer p : reversePrimers){
         double delta = Math.abs(targetTemp-p.getMt());
         if (delta > mtrange){
         //check that there's one better primer for this oligo
         for (Primer p2 : reversePrimers){
         if (!p2.equals(p) && p2.oligo.equals(p.oligo)){
         double delta2 = Math.abs(targetTemp-p2.getMt());
         if (delta > delta2){
         reversePrimers.remove(p);
         break;
         }
         }
         }
         }
         }*/
        ArrayList<ArrayList<Primer>> sets = populateReverseSets(acceptedReverse, pool);

        sets.add(0, acceptedForward);//put the forward primers in the first position
        return sets;
    }

    /**
     * Get a forward primer for the modified and unmodified version of each
     * target
     *
     * @param pool
     * @return
     */
    public ArrayList<Primer> generateForwardSet(List<Oligo> pool) {
        ArrayList<Primer> forwardPrimers = new ArrayList();
        for (Oligo oligo : pool) {
            //get all of the forward primers
            //5'->3' modified
            //5'->3' wildtype
            //3'->5' modified
            //3'->5' wildtype
            forwardPrimers.add(getModifiedForwardPrimer(oligo));
            forwardPrimers.add(getUnmodifiedForwardPrimer(oligo));
            forwardPrimers.add(getUnmodifiedAntisenseForwardPrimer(oligo));
            forwardPrimers.add(getModifiedAntisenseForwardPrimer(oligo));
        }
        //Forward primers have fewer options. Determine the ideal temp and
        //orientation of those first.
        //get the interquartile median of the melting temps (to ignore outliers)
        Melt.setMTs(forwardPrimers);
        if (!forceTemp) {
            Double[] mts = new Double[forwardPrimers.size()];
            for (int i = 0; i < forwardPrimers.size(); i++) {
                mts[i] = forwardPrimers.get(i).getMt();
            }
            targetTemp = iqmean(mts);
        }
        ArrayList<Primer> acceptedWT = new ArrayList<Primer>();
        ArrayList<Primer> acceptedMod = new ArrayList<Primer>();
        //get the best forward primer for each oligo
        for (Oligo oligo : pool) {
            Primer wtsense = null;
            Primer wtanti = null;
            Primer msense = null;
            Primer manti = null;
            for (Primer p : forwardPrimers) {
                if (p.oligo.equals(oligo)) {
                    if (p.sense && p.modified) {
                        msense = optimizePrimer(p);
                    } else if (p.sense && !p.modified) {
                        wtsense = optimizePrimer(p);
                    } else if (!p.sense && p.modified) {
                        manti = optimizePrimer(p);
                    } else if (!p.sense && !p.modified) {
                        wtanti = optimizePrimer(p);
                    }
                }
            }
            //choose the best unmodified forward primer
            if (wtanti == null && wtsense != null) {
                acceptedWT.add(wtsense);
            } else if (wtanti != null && wtsense == null) {
                acceptedWT.add(wtanti);
            } else if (Math.abs(targetTemp - wtsense.getMt()) <= Math.abs(targetTemp - wtanti.getMt())) {
                acceptedWT.add(wtsense);
            } else {
                acceptedWT.add(wtanti);
            }
            //choose the best mutant forward primer
            if (manti == null && msense != null) {
                acceptedMod.add(msense);
            } else if (manti != null && msense == null) {
                acceptedMod.add(manti);
            } else if (Math.abs(targetTemp - msense.getMt()) <= Math.abs(targetTemp - manti.getMt())) {
                acceptedMod.add(msense);
            } else {
                acceptedMod.add(manti);
            }
        }
        //correct the primers that are too far away from the target
        acceptedWT.addAll(acceptedMod);
        return acceptedWT;
    }

    //create a Primer that's one base longer
    public Primer extend(Primer p) {
        int newlen = p.seq.length() + 1;
        return changeLength(p, newlen);
    }

    //create a Primer that's one base shorter
    public Primer shorten(Primer p) {
        int newlen = p.seq.length() - 1;
        return changeLength(p, newlen);
    }

    private Primer changeLength(Primer p, int newlen) {
        if ((newlen > defaultPrimerLength + maxshift) || newlen < defaultPrimerLength - maxshift) {
            //don't change anything
            return p;
        } else {
            primerlength = newlen;
            if (p.forward) {//only valid to extend these from the end
                if (p.sense) {
                    if (p.modified) {
                        return getModifiedForwardPrimer(p.oligo);
                    } else {
                        return getUnmodifiedForwardPrimer(p.oligo);
                    }
                } else {
                    if (p.modified) {
                        return getModifiedAntisenseForwardPrimer(p.oligo);
                    } else {
                        return getUnmodifiedAntisenseForwardPrimer(p.oligo);
                    }
                }
            } else {
                Primer fromEnd = null;
                Primer fromStart = null;
                if (p.sense) {
                    fromEnd = getDownstreamReversePrimer(p.oligo, p.start, p.amplicon);
                    fromStart = getDownstreamReversePrimer(p.oligo, p.start - 1, p.amplicon);
                } else {
                    fromEnd = getUpstreamReversePrimer(p.oligo, p.start, p.amplicon);
                    fromStart = getUpstreamReversePrimer(p.oligo, p.start - 1, p.amplicon);
                }
                //compare the MTs, get the one closest to the target temp
                Melt.setMTs(new Primer[]{fromEnd,fromStart});
                
                Double dend = Math.abs(targetTemp -fromEnd.getMt());
                Double dstart = Math.abs(targetTemp - fromStart.getMt());
                if (dstart < dend) {
                    return fromStart;
                } else {
                    return fromEnd;
                }
            }
        }
    }

    private ArrayList<ArrayList<Primer>> populateReverseSets(ArrayList<Primer> reversePrimers,
            List<Oligo> pool) {
        ArrayList<Primer[]> sets = new ArrayList<>();
        for (int i = 0; i < Math.ceil((float) pool.size() / (float) ampliconLengths.size()); i++) {
            Primer[] set = new Primer[ampliconLengths.size()];
            sets.add(set);
        }
        //store wether an oligo is represented in the final reverse primer sets
        HashMap<Oligo, Boolean> isPresent = new HashMap<>();
        for (Oligo o : pool) {
            isPresent.put(o, false);
        }
        //use every other position in ampliconlengths first, so a non-full set will
        //increase resolution
        Integer[] ampPriority = new Integer[ampliconLengths.size()];
        for (int i = 0; i < ampPriority.length; i++) {
            ampPriority[i] = (int) Math.floor(i / 2);
            if (i % 2 == 1) {
                ampPriority[i] += (int) Math.floor((ampPriority.length + 1) / 2);
            }
        }
        List<Integer> ampOrder = Arrays.asList(ampPriority);

        //for each primer, if its oligo is not spoken for, crawl over primer sets 
        //seeing if it's a better fit for its position. If so, displace that 
        //primer. If it gets to the end of all sets, create a new set.
        int counter = 0;
        while ((isPresent.values().contains((Boolean) false)) && (counter < Math.pow(pool.size(), 2))) {
            System.out.println("DEBUG: " + counter);
            boolean unaltered = true; //no changes have been made during this loop
            //this will result in another set being added
            for (Integer priority = 0; priority < ampOrder.size(); priority++) {
                for (Primer p : reversePrimers) {
                    int ampIdx = ampliconLengths.indexOf(p.amplicon);//where does this primer go within its set?
                    int setIdx = ampOrder.get(ampIdx); //in which order are sets filled up?
                    if (setIdx == priority && (!isPresent.get(p.oligo))) {
                        //boolean placed = false;
                        for (Primer[] arr : sets) {
                            if (arr[ampIdx] != null) { //something is already here
                                Primer p2 = arr[ampIdx];
                                double delta = Math.abs(targetTemp - p.getMt());
                                double delta2 = Math.abs(targetTemp - p2.getMt());
                                if (delta < delta2) { //replace it
                                    isPresent.put(p2.oligo, false);
                                    //placed = true;
                                    isPresent.put(p.oligo, true);
                                    arr[ampIdx] = p;
                                    unaltered = false;
                                    break;
                                }
                            } else { //the space is empty
                                //placed = true;
                                isPresent.put(p.oligo, true);
                                arr[ampIdx] = p;
                                unaltered = false;
                                break;
                            }
                        }
                    } //cover the case where the oligo does have a primer present
                    //but a better primer fits elsewhere without displacement
                    else if (setIdx == priority && (isPresent.get(p.oligo))) {
                        //get the other primer
                        Primer other = null;
                        Integer otherSetIdx = null;
                        outerloop:
                        for (int i = 0; i < sets.size(); i++) {
                            Primer[] arr = sets.get(i);
                            for (Primer p2 : arr) {
                                if (p2 != null && p2.amplicon != p.amplicon && p.oligo.equals(p2.oligo)) {
                                    other = p2;
                                    otherSetIdx = i;
                                    break outerloop;
                                }
                            }
                        }
                        //does the new primer have a better temp?
                        if (other != null && Math.abs(targetTemp - p.getMt()) < Math.abs(targetTemp - other.getMt())) {
                            //if there's an unoccupied slot, replace
                            for (Primer[] arr : sets) {
                                if (arr[ampIdx] == null) {
                                    arr[ampIdx] = p;
                                    int otherIdx = ampliconLengths.indexOf(other.amplicon);
                                    Primer[] adjusted = sets.get(otherSetIdx);
                                    adjusted[otherIdx] = null;
                                    sets.set(otherSetIdx, adjusted);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            //if there is some oligo that can't fit
            if (unaltered && isPresent.values().contains((Boolean) false)) {
                Primer[] set = new Primer[ampliconLengths.size()];
                sets.add(set);
            }
            counter++;
        }
        //DEBUG
        /*for (int i = 0; i < sets.size(); i++){
         Primer[] arr = sets.get(i);
         String s = "";
         for (Primer p: arr){
         //System.out.println(p);
         if (p != null){
         s = s.concat(p.oligo.name + "." + p.amplicon + " " );
         }
         else{
         s = s.concat("null ");
         }
         }
         System.out.println(s);
         }*/

        //convert to ArrayLists
        ArrayList<ArrayList<Primer>> arr = new ArrayList();
        for (Primer[] set : sets) {
            ArrayList<Primer> a = new ArrayList(Arrays.asList(set));
            //wipe out the nulls
            List nl = new ArrayList();
            nl.add(null);
            a.removeAll(nl);

            arr.add(a);
        }
        return arr;
    }

    //extract the sequences from each primer and return their melting temps
    //in the same order
    public static Double[] getMeltingTemps(ArrayList<Primer> primers) {
        List<String> seqs = new ArrayList();
        for (Primer p : primers) {
            seqs.add(p.seq);
        }
        Double[] res = new Double[0];
        try {
            res = Melt.getMT(seqs);
        } catch (IOException ex) {
            Logger.getLogger(PCR.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(PCR.class.getName()).log(Level.SEVERE, null, ex);
        }
        return res;
    }

    /**
     * Alter the primers in the list so that for each oligo, all falling outside
     * the optimal range are filtered out. Unless it is not possible to create a
     * primer in the correct range, then get as close as you can
     *
     * @param list
     * @return
     */
    public ArrayList<Primer> correctOutlyingPrimers(ArrayList<Primer> list) {
        return correctOutlyingPrimers(list, forceTemp);
    }

    //give the ability to temporarily override forceTemp, as is done when optimizing
    //a single primer
    private ArrayList<Primer> correctOutlyingPrimers(ArrayList<Primer> list, boolean lockTemp) {
        int counter = 0;
        boolean running = true;
        while (running && counter < maxshift) {
            counter++;
            list = evaluateAndShiftPrimers(list, lockTemp);
        }
        return list;
    }

    /**
     * Alter the length of the given primer to maximize it relative to the
     * current target temp
     *
     * @param p
     * @return
     */
    public Primer optimizePrimer(Primer p) {
        ArrayList<Primer> arr = new ArrayList(1);
        arr.add(p);
        ArrayList<Primer> res = correctOutlyingPrimers(arr, true);
        return res.get(0);
    }

    /**
     * If an oligo has any primers in the acceptable temperature range, retain
     * only those. If it does not, shift each primer's length by 1. This is
     * called multiple times by correctOutlyingPrimers
     *
     * @param list
     * @return
     */
    public ArrayList<Primer> evaluateAndShiftPrimers(ArrayList<Primer> list, boolean lockTemp) {
        if (!lockTemp) {
            Double[] mt = getMeltingTemps(list);
            //get the weighted average temperature
            targetTemp = iqmean(mt);
        }

        //get all Oligos
        ArrayList<Oligo> oligos = new ArrayList();
        for (Primer p : list) {
            Oligo o = p.oligo;
            if (!oligos.contains(o)) {
                oligos.add(o);
            }
        }

        //if any oligo doesn't have a primer within an acceptable range, attempt
        //to create a new set of primers by varying the length
        ArrayList<Primer> allAccepted = new ArrayList();
        for (Oligo o : oligos) {
            ArrayList<Primer> primers = new ArrayList();
            ArrayList<Primer> accepted = new ArrayList();
            for (Primer p : list) {
                if (o.equals(p.oligo)) {
                    primers.add(p);
                    //System.out.println("DEBUG: mt = " + p.getMt());
                    if (Math.abs(p.getMt() - targetTemp) <= mtrange) {
                        accepted.add(p);
                    }
                }
            }
            if (accepted.isEmpty()) {
                for (int i = 0; i < primers.size(); i++) {
                    Primer p = primers.get(i);
                    //compare to the effect if one base longer or shorter
                    double delta = Math.abs(targetTemp - p.getMt());
                    Primer shorter = shorten(p);
                    double shortDelta = Math.abs(targetTemp - shorter.getMt());
                    Primer longer = extend(p);
                    double longDelta = Math.abs(targetTemp - longer.getMt());

                                //DEBUG
                    /*ArrayList<Primer> debugging = new ArrayList(3);
                     debugging.add(p);
                     debugging.add(shorter);
                     debugging.add(longer);
                     for (Primer d : debugging){
                     String s = "";
                     s = s.concat(d.oligo.name + "\t");
                     s = s.concat("Fwd: " + d.forward + "\t");
                     s = s.concat("Sense: " + d.sense + "\t");
                     s = s.concat("Mod: " + d.modified + "\t");
                     s = s.concat("Mt: " + d.getMt() + "\t");
                     s = s.concat("Len: " + d.seq.length() +"\t");
                     s = s.concat(d.seq);
                     System.out.println(s);
                     }*/
                    //pick the alternate that did best
                    Primer candidate = longer;
                    if (Math.min(shortDelta, longDelta) == shortDelta) {
                        candidate = shorter;
                    }
                    //if the original did better
                    if (delta <= Math.min(shortDelta, longDelta)) {
                        accepted.add(p); //not going to get any better
                    } //if the candidate did better
                    else {
                        accepted.add(candidate);
                    }
                }
            }
            allAccepted.addAll(accepted);
        }
        return allAccepted;
    }

    /**
     * determines the ideal temperature for a pool of primers using one primer
     * per oligo, determined by least means squared
     *
     * @param list
     * @return
     */
    //public static double lmsTemp(ArrayList<Primer> list){
    //  Double[] mts = getMeltingTemps(list);
    //}
    //convenience method to make a forward primer
    public static Primer forwardPrimer(String seq, Oligo oligo, int start, boolean sense, boolean modified) {
        Primer p = new Primer(seq, oligo, 0, start, true, sense, modified);
        
        return p;
    }

    /**
     * Calculate the interquartile mean of the given list as a quick and dirt
     * way of ignoring outliers. If n<4, simply return the arithmatic mean
     *
     * @param list
     * @return
     */
    private static Double iqmean(Double[] list) {
        int n = list.length;
        if (n < 4) {
            return mean(list);
        } else {
            int lowidx = ((Double) Math.floor(n / 4)).intValue();
            int highidx = ((Double) Math.ceil(n * 3 / 4)).intValue();
            Double[] sublist = new Double[highidx - lowidx];
            for (int i = lowidx; i < highidx; i++) {
                sublist[i - lowidx] = list[i];
            }
            return mean(sublist);
        }
    }

    private static Double mean(Double[] list) {
        Double sum = 0.0;
        for (Double d : list) {
            sum += d;
        }
        return Double.valueOf(sum / list.length);
    }
}
