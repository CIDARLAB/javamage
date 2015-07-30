package mage.Tools.Pcr;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;

/**
 * Class to generate MASC PCR primers for oligos see
 * http://wanglab.c2b2.columbia.edu/publications/2011_MIE_Wang.pdf
 *
 * Assumptions Begin the primer placing the first modified base at the 3'
 * position Direction doesn't matter- if the opposite strand gets better score,
 * use that Genome input is already 5' to 3' Uses a complete, circular genome
 * Primer for deletion: put the start 1/4 * primerLength before the deletion
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
 * @author mquintin
 *
 */
public class PCR {
    public String genome;
    
    private PrimerFactory pf;

    private List<Integer> ampliconLengths = Arrays.asList(100, 150, 200, 250, 300, 400, 500, 700, 850);
    //legal range of primer lengths
    private static double mtrange = 1.5;//how far from the final temp is OK?

    public static boolean forceTemp = true; //if true, primers are based on the targetTemp
    //if false, primers are based on the average of the forward primer temps
    public static double targetTemp = 60.0;//default target melting temp in degrees C
    private int primerLength;

    private double shiftRange;
    private int lenMin;
    private int lenMax;
    private static final double defaultShiftRange = 0.15;
    private static final int defaultLenMin = 16;
    private static final int defaultLenMax = 30;
    

    public PCR() throws IOException{
        this(defaultShiftRange, defaultLenMin, defaultLenMax);
    }
    
    public PCR(String genome){
        this(genome, defaultShiftRange, defaultLenMin, defaultLenMax, null);
    }
    
    //defaults to loading the genome file in the working directory
    public PCR(double shiftRange, int lenMin, int lenMax) throws IOException{
        this(FASTA.readFFN(Oligo.Directory, Oligo.Genome), shiftRange, lenMin, lenMax, null);
        System.err.println("Merlin-DEBUG: defaulting to oligo file " + Oligo.Directory + "/" + Oligo.Genome);
    }
    
    public PCR(String genome, double shiftRange, int lenMin, int lenMax, List<Integer> ampliconLengths){
        if (ampliconLengths != null) this.ampliconLengths = ampliconLengths;
        if (Melt.directory == null) {
            Melt.directory = Oligo.Directory;
        }
        this.genome = genome;
        this.pf = new PrimerFactory(genome);
        this.lenMin = lenMin;
        this.lenMax = lenMax;
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
                    reversePrimers.add(pf.getDownstreamReversePrimer(o, amp));
                } else {
                    reversePrimers.add(pf.getUpstreamReversePrimer(o, amp));
                }
            }
        }
        ArrayList<Primer> acceptedReverse = optimizePrimers(reversePrimers);
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
    public ArrayList<Primer> generateForwardSet(List<Oligo> pool) throws IOException {
        ArrayList<Primer> forwardPrimers = new ArrayList();
        for (Oligo oligo : pool) {
            //get all of the forward primers
            //5'->3' modified
            //5'->3' wildtype
            //3'->5' modified
            //3'->5' wildtype
            forwardPrimers.add(pf.getModifiedForwardPrimer(oligo));
            forwardPrimers.add(pf.getUnmodifiedForwardPrimer(oligo));
            forwardPrimers.add(pf.getUnmodifiedAntisenseForwardPrimer(oligo));
            forwardPrimers.add(pf.getModifiedAntisenseForwardPrimer(oligo));
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
            //choose the best mutant forward primer and its corresponding unmodified primer
            if (manti == null && msense != null) {
                acceptedMod.add(msense);
                acceptedWT.add(wtsense);
            } else if (manti != null && msense == null) {
                acceptedMod.add(manti);
                acceptedWT.add(wtanti);
            } else if (Math.abs(targetTemp - msense.getMt()) <= Math.abs(targetTemp - manti.getMt())) {
                acceptedMod.add(msense);
                acceptedWT.add(wtsense);
            } else {
                acceptedMod.add(manti);
                acceptedWT.add(wtanti);
            }
        }
        //correct the primers that are too far away from the target
        acceptedWT.addAll(acceptedMod);
        return acceptedWT;
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
            //System.out.println("DEBUG: " + counter);
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
    
    public ArrayList<Primer> optimizePrimers(ArrayList<Primer> list){
        ArrayList<Primer> newList = new ArrayList();
        for (Primer p : list){
            newList.add(optimizePrimer(p));
        }
        return newList;
    }

    /**
     * Alter the length of the given primer to maximize it relative to the
     * current target temp
     *
     * @param p
     * @return
     */
    public Primer optimizePrimer(Primer p) {
        Primer opt = null;
        try {
            Anneal anneal = new Anneal(genome, targetTemp, p.oligo, p.amplicon,
                    shiftRange, lenMin, lenMax, p.forward, p.sense, p.modified);
            opt = anneal.getOptimizedPrimer();
        } catch (IOException e) {
            opt = p;
        }
        return opt;
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
    
    //Getters and Setters

    public static double getMtrange() {
        return mtrange;
    }

    public static void setMtrange(double mtrange) {
        PCR.mtrange = mtrange;
    }

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
    
    public String getGenome() {
        return genome;
    }

    public void setGenome(String genome) {
        this.genome = genome;
    }

    public PrimerFactory getPf() {
        return pf;
    }

    public void setPf(PrimerFactory pf) {
        this.pf = pf;
    }

    public int getPrimerLength() {
        return primerLength;
    }

    public void setPrimerLength(int primerLength) {
        this.primerLength = primerLength;
    }

    public double getShiftRange() {
        return shiftRange;
    }

    public void setShiftRange(double shiftRange) {
        this.shiftRange = shiftRange;
    }

    public int getLenMin() {
        return lenMin;
    }

    public void setLenMin(int lenMin) {
        this.lenMin = lenMin;
    }

    public int getLenMax() {
        return lenMax;
    }

    public void setLenMax(int lenMax) {
        this.lenMax = lenMax;
    }
    
    public List<Integer> getAmpliconLengths() {
        return ampliconLengths;
    }

}
