/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import mage.Core.Oligo;
import mage.Tools.FASTA;
import mage.Tools.OligoStats;
import org.biojava3.core.sequence.DNASequence;
import static test.TestReplacementPrediction.oligoID;

/**
 *
 * @author mquintin
 */
public class TestDiversityPrediction {
    static ArrayList<Oligo> oligos = new ArrayList(0);
    static String genomeDir = "C:/Users/mquintin/OneDrive/github/javamage/testing/";
    static String genomeFile = "NC000913.fasta";
    static String[] names = new String[0];
    static Boolean[] sense = new Boolean[0];
    static Integer[] replicore = new Integer[0];
    static Integer[] lp = new Integer[0];
    static Integer[] rp = new Integer[0];
    static String[] mut = new String[0];
    static String[] seq = new String[0];
    static ArrayList<Integer> ignore = new ArrayList(0);
    
    public static void populateOligos() throws Exception{
        names = new String[]{"ompT","lon","endA","ilVG-f","argR","iclRp1"};
        sense = new Boolean[]{false,true,true,true,true,false};
        replicore = new Integer[]{1,1,2,1,2,1};
        lp = new Integer[]{585624,458906,3090358,3951536,3384911,4223694};
        rp = new Integer[]{585624,458906,3090358,3951536,3384911,4223694};
        mut = new String[]{"M","M","M","I","M","M"};
        seq = new String[]{"T","T","A","AT","T","G"};
        
        //ignore.add(4);
        //ignore.add(1);
        
        String genome = FASTA.readFFN(genomeDir, genomeFile);
        
        int i = 0;
        while (i < names.length){
            //System.out.println(seq[i]+" "+lp[i]+" "+rp[i]+" "+replicore[i]+" "+sense[i]+" "+names[i]);
            if (!ignore.contains(i)){
                if ("M".equals(mut[i])){
                    //System.out.println("Creating " + names[i]);
                    oligos.add(Oligo.MismatchFactory(genome, seq[i], 
                            lp[i], rp[i], replicore[i], sense[i], names[i]));
                }
                else{
                    //System.out.println("Creating " + names[i]);
                    oligos.add(Oligo.InsertionFactory(genome, seq[i], 
                            lp[i], replicore[i], sense[i], names[i]));
                }
            }
            i++;
        }
        
        //System.out.println(oligos);
        System.out.println("getAggregateAnyARE: " + OligoStats.getAggregateAnyARE(oligos));
	System.out.println("getAggregateSumARE: " + OligoStats.getAggregateSumARE(oligos));
        
        double prob = OligoStats.getAggregateAnyARE(oligos);
        List<Integer> loci = Arrays.asList(new Integer[]{0,1,2,3,4,5,6});
        int m = Collections.max(loci);
        for (int c : loci){
            double res = OligoStats.lociProbability(m, c, prob, 1);
            System.out.println(c + " out of " + m + " loci: " + res);
        }
        Object res = OligoStats.getDiscreteDiversityTable(oligos, 1);
		System.out.println("Diversity table 1:");
		System.out.println(res);
		System.out.println();
                System.out.println();
                System.out.println(OligoStats.probabilitiesAfterCycles(oligos, 1));

    }
    
    public static void main(String[] args) throws Exception{
        populateOligos();
        //testBioJava();
    }
    
    public static void testBioJava() throws Exception{
        String seq = "TAGAGGCGCTACCGAACGGTGCCGCAGAAACCGCGCCCTGACGGGTTAACATGCCTTCGATTT"+
"GCACCACGCTAAGACCCGGTACAAACGGTGCCAAATGCGCTTTCACGCCGATCGGTCCCATACCCGGACC"+
"ACCACCGCCGTGCGG";
        String pre = "TAGAGGCGCTACCGAACGGTGCCGCAGAAACCGCGCCCTGACGGGTTAACATGCCTTCGATTTGCACCACGCTA";
        String post = "GACCCGGTACAAACGGTGCCAAATGCGCTTTCACGCCGATCGGTCCCATACCCGGACCACCACCGCCGTGCGG";
        String target = "A";
        DNASequence dna = new DNASequence(pre+target+post);
        System.out.println(dna);
        String genome = FASTA.readFFN(genomeDir, genomeFile);
        Oligo o = Oligo.MismatchFactory(genome,target,2090358,3090358,2,true,"testEndA");
        System.out.println(mage.Tools.SequenceTools.ReverseCompliment(o.sequence));
    }
}
