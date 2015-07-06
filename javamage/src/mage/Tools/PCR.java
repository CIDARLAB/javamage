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
 * TODO: Confirm method for forward unmodified deletion primer: is currently only including 5 bases into the cut
 * This should be redone to put the leftmost deleted base at the end
 * 
 * TODO: Confirm method for modified antisense forward primer. If the 3' buffer length 
 * is shorter than the primer length, primer may run off the oligo span. Need to
 * add part of the genome seq to cover this case.
 * 
 * @author mquintin
 *
 */

//TODO: Account for cases where the primers may extend past the start/end of the genome 
public class PCR {
	
    private String genome;
    private int primerlength;
    private static List<Integer> ampliconLengths = Arrays.asList(100,150,200,250,300,400,500,700,850);
    private static int defaultPrimerLength = 20; //generally 18-30
    //number of valid primers to explore, keeping the first modified base in the 3' half
    private static double maxshift = Math.floor(defaultPrimerLength/2.0);
    private static double mtrange = 0.5;//how far from the final temp is OK?

    public PCR(String genome){
        this.genome = genome;
        this.primerlength = defaultPrimerLength;
    }
    
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
    public Primer getUnmodifiedForwardPrimer(Oligo oligo){
            //position of the last targeted base on the reference genome
            int end = oligo.getGenomeStart() + oligo.target_position + oligo.target_length;

            //for deletions, try to cover as much of the deleted region as possible
            if (oligo.getOligoType() == OligoType.DELETION){
                    end = oligo.getGenomeEnd() + oligo.target_position - (oligo.getLength() + oligo.getGenomeStart());
                    end = oligo.target_position + oligo.getGenomeEnd() - oligo.getLength();
            }
            end = end - 1; //necessary to align this with the modified primers 
            String seq = genome.substring(end - primerlength, end);
            Primer primer = forwardPrimer(seq, oligo, end - primerlength, true, false);
            return primer;
    }

    public Primer getUnmodifiedAntisenseForwardPrimer(Oligo oligo){
            //position of the last targeted base on the reference genome
            int start = oligo.getGenomeStart() + oligo.target_position;

            String seq = genome.substring(start, primerlength + start);
            seq = SequenceTools.ReverseCompliment(seq);
            Primer primer = forwardPrimer(seq, oligo, start, false, false);
            return primer;
    }

    /**Get the forward PCR primer for the modified sequence, given the Oligo
     * 
     * @param genome
     * @param start
     * @return
     */
    public Primer getModifiedForwardPrimer(Oligo oligo){
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
            Primer primer = forwardPrimer(seq, oligo, start, true, true);
            return primer;
    }

    /**get the primer in the opposite direction. In this case, the last modified
    *base in the oligo will be the first base in the primer
    */
    public Primer getModifiedAntisenseForwardPrimer(Oligo oligo){
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
            seq = SequenceTools.ReverseCompliment(seq);
            Primer primer = forwardPrimer(seq, oligo, start, true, true);
            return primer;
    }


    /**Get a single PCR reverse primer
     * 
     * @param oligo
     * @param start start position for the forward primer on the reference genome
     * @param amplength length of the amplicon, NOT the primer
     * @return primer sequence as a String of length = primerlength
     */
    public Primer getDownstreamReversePrimer(Oligo oligo, int start, int amplength){
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

            Primer primer = new Primer(seq, oligo, amplength, start, false, true);

            return primer;
    }

    public Primer getUpstreamReversePrimer(Oligo oligo, int start, int amplength){
            String gen = genome;
            //gen = new StringBuffer(gen).reverse().toString();
            //int s = gen.length()-start; //where to start this primer
            int s = start + oligo.target_length - amplength;
            //if the primer is too close to the end, it will run off
            //take enough bases from the start, paste them onto the end
            if (gen.length() - s < amplength){
                    gen = gen.concat(gen.substring(0,amplength));
            }
            if (s < 0){ //handle spanning the origin
                gen = gen.substring(gen.length()+s) + gen.substring(0,amplength);
                s = 0;
            }
            String seq = gen.substring(s, s + primerlength);
            //seq = SequenceTools.ReverseCompliment(seq);

            Primer primer = new Primer(seq, oligo, amplength, start, false, true);

            return primer;
    }


    /**Get the unmodified forward (first element), modified forward (second),
     *  and reverse (third through eleventh elements) 
     * PCR primers for this oligo
     * 
     * @param oligo
     * @param genome
     * @return
     */
    public List<String> getMASCPCRPrimers(Oligo oligo,String genome){
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
    public Primer getSenseReversePrimer(Oligo oligo, int len, int baselength){
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

        Primer primer = new Primer(seq, oligo, baselength, start, false, true);

        return primer;
    }

    public Primer getAntisenseReversePrimer(Oligo oligo, int len, int baselength){
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

        Primer primer = new Primer(seq, oligo, baselength, start, false, false);

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
    public List<List<String>> getMASCPCRPrimersV0(List<Oligo> pool) throws IOException{
            // First we read in the genome
            String genome = FASTA.readFFN(Oligo.Directory,Oligo.Genome);

            List<List<String>> primerset = new ArrayList<List<String>>();
            for (Oligo oligo:pool){
                    primerset.add(getMASCPCRPrimers(oligo,genome));
            }
            return primerset;
    }

    /**get a list of lists of primers. The first list is the forward primers, each
     * subsequent list is a set of reverse primers that is grouped together for
     * one MASC PCR reaction. The forward primers should be paired with these sets
     * by checking the Primer.oligo field
     * 
     * @param pool a List of Oligos
     * @return
     * @throws IOException 
     */
    public ArrayList<ArrayList<Primer>> generateAllPrimers(List<Oligo> pool) throws IOException{
        ArrayList<Primer> forwardPrimers = new ArrayList();
        ArrayList<Primer> reversePrimers = new ArrayList();
        double range = 0.05; //how much variablility is acceptible in the position?
            //ex: if range = .05 the amplicon must be from 95% to 105% of the given length

        for (Oligo oligo : pool){
            //get the forward primers
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
        Double[] mts = Melt.getMT(forwardPrimers);
        Double targetTemp = iqmean(mts);
        ArrayList<Primer> acceptedForward = new ArrayList<Primer>();
        //get the best forward primer for each oligo
        for (Oligo oligo : pool){
            Primer wtsense = null;
            Primer wtanti = null;
            Primer msense = null;
            Primer manti = null;
            for (Primer p : forwardPrimers){
                if (p.oligo.equals(oligo)){
                    if (p.sense && p.modified){
                        msense = p;
                    }
                    else if (p.sense && !p.modified){
                        wtsense = p;
                    }
                    else if (!p.sense && p.modified){
                        manti = p;
                    }
                    else if (!p.sense && !p.modified){
                        wtanti = p;
                    }
                }
            }
            //choose the best unmodified forward primer
            if (wtanti == null && wtsense != null){
                acceptedForward.add(wtsense);
            }
            else if (wtanti != null && wtsense == null){
                acceptedForward.add(wtanti);
            }
            else if (Math.abs(targetTemp-wtsense.getMt()) <= Math.abs(targetTemp-wtanti.getMt())){
                acceptedForward.add(wtsense);                
            }
            else{
                acceptedForward.add(wtanti);
            }
            //choose the best mutant forward primer
            if (manti == null && msense != null){
                acceptedForward.add(msense);
            }
            else if (manti != null && msense == null){
                acceptedForward.add(manti);
            }
            else if (Math.abs(targetTemp-msense.getMt()) <= Math.abs(targetTemp-manti.getMt())){
                acceptedForward.add(msense);                
            }
            else{
                acceptedForward.add(manti);
            }
        }
        //correct the primers that are too far away from the target
        acceptedForward = correctOutlyingPrimers(acceptedForward);
        double adjustedTargetTemp = iqmean(Melt.getMT(acceptedForward));
        int counter = 0;
        while (adjustedTargetTemp != targetTemp && counter < maxshift){
            counter++;
            targetTemp = adjustedTargetTemp;
            acceptedForward = correctOutlyingPrimers(acceptedForward,targetTemp);
            adjustedTargetTemp = iqmean(Melt.getMT(acceptedForward));
        }

        //get the reverse primers that correspond with the selected forwards
        for (int amp : ampliconLengths){
            for (Primer p : acceptedForward){
                Oligo o = p.oligo;
                if (p.sense){
                    reversePrimers.add(getDownstreamReversePrimer(o, p.start, amp));
                }
                else{
                    reversePrimers.add(getUpstreamReversePrimer(o, p.start, amp));
                }
            }       
        }
        ArrayList<Primer> acceptedReverse = correctOutlyingPrimers(reversePrimers,targetTemp);
        counter = 0;
        while (counter < maxshift){
            if(reversePrimers.equals(acceptedReverse)){
                break;
            }
            reversePrimers = acceptedReverse;
            acceptedReverse = correctOutlyingPrimers(reversePrimers,targetTemp);
            counter++;
        }
        //remove any reverse primer that is still unacceptable
        for(Primer p : reversePrimers){
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
        }
        ArrayList<ArrayList<Primer>> sets = populateReverseSets(reversePrimers,pool,targetTemp);
        
        sets.add(0, acceptedForward);//put the forward primers in the first position
        return sets;
    }
    
    private ArrayList<Primer> correctOutlyingPrimers(ArrayList<Primer> primers, double targetTemp){
        //TODO
        return primers;
    }
    
    //create a Primer that's one base longer
    public Primer extend(Primer p, double targetTemp){
            int newlen = p.seq.length() +1;
        return changeLength(p, targetTemp, newlen);  
    }
    
    //create a Primer that's one base shorter
    public Primer shorten(Primer p, double targetTemp){
        int newlen = p.seq.length() -1;
        return changeLength(p, targetTemp, newlen);        
    }
    
    private Primer changeLength(Primer p, double targetTemp, int newlen){
        if ((newlen > defaultPrimerLength + maxshift) || newlen < defaultPrimerLength - maxshift){
            //don't change anything
            return p;
        }
        else{
            primerlength = newlen;
            if (p.forward){//only valid to extend these from the end
                if (p.sense){
                    if (p.modified){
                        return getModifiedForwardPrimer(p.oligo);
                    }
                    else{
                        return getUnmodifiedForwardPrimer(p.oligo);
                    }
                }
                else{
                    if (p.modified){
                        return getModifiedAntisenseForwardPrimer(p.oligo);
                    }
                    else{
                        return getUnmodifiedAntisenseForwardPrimer(p.oligo);
                    }
                }
            }
            else{
                Primer fromEnd = null;
                Primer fromStart = null;
                if (p.sense){
                    fromEnd = getDownstreamReversePrimer(p.oligo, p.start, p.amplicon);
                    fromStart = getDownstreamReversePrimer(p.oligo, p.start-1, p.amplicon);
                }
                else{
                    fromEnd = getUpstreamReversePrimer(p.oligo, p.start, p.amplicon);
                    fromStart = getUpstreamReversePrimer(p.oligo, p.start-1, p.amplicon);
                }
                //compare the MTs, get the one closest to the target temp
                String[] seqs = new String[2];
                seqs[0] = fromEnd.seq;
                seqs[1] = fromStart.seq;
                Double[] mts = new Double[2];
                try {
                    mts = Melt.getMT(seqs);
                } catch (IOException | InterruptedException ex) {
                    Logger.getLogger(PCR.class.getName()).log(Level.SEVERE, null, ex);
                    mts[0] = 0.0;
                    mts[1] = 0.0;
                }
                Double dend = Math.abs(targetTemp - mts[0]);
                Double dstart = Math.abs(targetTemp - mts[1]);
                if (dstart < dend){
                    return fromStart;
                }
                else{
                    return fromEnd;
                }
            }
        }
    }
    
    private ArrayList<ArrayList<Primer>> populateReverseSets(ArrayList<Primer> reversePrimers, 
            List<Oligo> pool,  double targetTemp){
        ArrayList<ArrayList<Primer>> sets = new ArrayList();
        for (int i = 0; i < Math.ceil(pool.size()/ampliconLengths.size()); i++){
            sets.add(new ArrayList());
        }
        //store wether an oligo is represented in the final reverse primer sets
        HashMap<Oligo,Boolean> isPresent = new HashMap<>();
        for (Oligo o : pool){
            isPresent.put(o, false);
        }
        //use every other position in ampliconlengths first, so a non-full set will
        //increase resolution
        int[] ampPriority = new int[ampliconLengths.size()];
        for (int i=0; i <= ampPriority.length/2; i += 2){
            ampPriority[i] = 2*i;
        }
        for (int i=1; i < ampPriority.length; i += 2){
            ampPriority[i] = (int) (Math.floor(ampPriority.length/2) + i); 
        }
        //for each primer, if its oligo is not spoken for, crawl over primer sets 
        //seeing if it's a better fit for its position. If so, displace that 
        //primer. If it gets to the end of all sets, create a new set.
        int counter = 0;
        while((isPresent.values().contains((Boolean) false))&&(counter < Math.pow(pool.size(),2))){
            for (int priority = 0; priority < ampPriority.length; priority++){
                for (Primer p : reversePrimers){
                    int ampIdx = ampliconLengths.indexOf(p.amplicon);
                    int setIdx = ampPriority[ampIdx];
                    if (ampPriority[setIdx]==priority && (!isPresent.get(p.oligo))){
                        boolean placed = false;
                        for (ArrayList<Primer> arr : sets){
                            if (arr.get(setIdx) != null){ //something is already here
                                Primer p2 = arr.get(setIdx);
                                double delta = Math.abs(targetTemp-p.getMt());
                                double delta2 = Math.abs(targetTemp-p2.getMt());
                                if (delta < delta2){ //replace it
                                    isPresent.put(p2.oligo,false);
                                    placed = true;
                                    isPresent.put(p.oligo,true);
                                    arr.set(setIdx,p);
                                }
                            } else { //the space is empty
                                placed = true;
                                isPresent.put(p.oligo,true);
                                arr.set(setIdx,p);
                            }
                        }
                        if (!placed){
                            ArrayList<Primer> newArr = new ArrayList<Primer>(ampliconLengths.size());
                            sets.add(newArr);
                        }
                    }
                }
            }
            counter++;
        }
        return sets;
    }

    //extract the sequences from each primer and return their melting temps
    //in the same order
    public static Double[] getMeltingTemps(ArrayList<Primer> primers){
        List<String> seqs = new ArrayList();
        for (Primer p : primers){
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


    public static ArrayList<Primer> correctOutlyingPrimers(ArrayList<Primer> list){
        Double[] mt = getMeltingTemps(list);
        //get the weighted average temperature

        //if any oligo doesn't have a primer within an acceptable distance, attempt
        //to create a new set of primers by varying the length

        //if a set could be created, replace the old set and rerun


        return list;
    }

    /**determines the ideal temperature for a pool of primers using one primer
     * per oligo, determined by least means squared
     * 
     * @param list
     * @return 
     */
    //public static double lmsTemp(ArrayList<Primer> list){
      //  Double[] mts = getMeltingTemps(list);
    //}


            //convenience method to make a forward primer
    public static Primer forwardPrimer(String seq, Oligo oligo, int start, boolean sense, boolean modified){
        return new Primer(seq,oligo,0,start,true,sense,modified);
    }
    
    /**Calculate the interquartile mean of the given list as a quick and dirt way
     * of ignoring outliers.
     * If n<4, simply return the arithmatic mean
     * 
     * @param list
     * @return 
     */
    private Double iqmean(Double[] list){
        int n = list.length;
        if ( n < 4){
            return mean(list);
        }
        else{
            int lowidx = ((Double) Math.floor(n/4)).intValue();
            int highidx = ((Double) Math.ceil(n*3/4)).intValue();
            Double[] sublist = new Double[highidx-lowidx];
            for (int i = lowidx; i < highidx; i++){
                sublist[i-lowidx] = list[i];
            }
            return mean(sublist);
        }
    }
    private Double mean(Double[] list){
        Double sum = 0.0;
        for (Double d : list){
            sum += d;
        }
        return Double.valueOf(sum/list.length);
    }
}
