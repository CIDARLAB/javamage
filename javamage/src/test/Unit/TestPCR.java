package test.Unit;

import java.io.IOException;
import java.util.ArrayList;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;
import mage.Tools.Pcr.PCR;
import test.Constants;

public class TestPCR {

	//Isaacs et al supplement put the modified base at the end (3') of the forward primer.
    //Should fix forward primer design to put modified bases towards the back of the sequence, not the start
    //confirmed from same paper that the forward primer is not reversed or complementary
    //here's a chunk of sequence we use from the E coli genome
    //remember substrings are built with 0 at the first index
    //80 TTCTGAACTG
    //90 GTTACCTGCC
    //100 GTGAGTAAAT
    //110 TAAAATTTTA
    public static void main(String[] args) throws Exception {
        //testGetUnmodifiedForwardPrimerFromOligo();
        //testGetModifiedForwardPrimerFromOligo();
        //testGetMASCPCRPrimersForOligo();
        //testPCRbyReplicatingResults();
        //testGeneratePrimerSets();
        //testGetForwardPrimers();
        testOptimizePrimer();
    }


    public static void testOptimizePrimer() throws IOException, Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        PCR pcr = new PCR(genome);
        Oligo o = Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2, true, "oligo10");
        Primer p = pcr.getPf().getModifiedForwardPrimer(o);
        Primer opt = pcr.optimizePrimer(p);
        System.out.println(opt.seq + " MT: " + opt.getMt());
                
    }

    public static void testGetForwardPrimers() throws IOException, Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        PCR pcr = new PCR(genome);
        ArrayList<Oligo> pool = new ArrayList<>();
        pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2, true, "oligo1"));
        pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2, true, "oligo2"));
        pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2, true, "oligo3"));
        pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2, true, "oligo4"));
        pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2, true, "oligo5"));
        pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2, true, "oligo6"));
        pool.add(Oligo.InsertionFactory(genome, "GACGATATGAAC", 7900, 2, true, "oligo7"));
        pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2, true, "oligo8"));
        pool.add(Oligo.InsertionFactory(genome, "ATTTGACCT", 15426, 2, true, "oligo9"));
        pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2, true, "oligo10"));

        ArrayList<Primer> arr = pcr.generateForwardSet(pool);
        //System.out.println(res);

        for (Primer p : arr) {
            String s = "";
            //System.out.println(p);
            if (p != null) {
                s = s.concat(p.oligo.name + "\t");
                s = s.concat("Fwd: " + p.forward + "\t");
                s = s.concat("Sense: " + p.sense + "\t");
                s = s.concat("Mod: " + p.modified + "\t");
                s = s.concat("Mt: " + p.mt + "\t");
                s = s.concat("Len: " + p.seq.length() + "\t");
                s = s.concat(p.seq);
                System.out.println(s);
            } else {
                System.out.println(p);
            }
        }
    }

    public static void testGetUnmodifiedForwardPrimerFromOligo() throws Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2, true, "in1");
        String expect = "GGTTACCTGCCGTGAGTAAA";
        //String iprimer = PCR.getUnmodifiedForwardPrimer(insert, genome);
        //System.out.println("Forward from insert oligo: " + iprimer);
        //System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.equals(iexpect)));
        PCR pcr = new PCR(genome);
        Primer p = pcr.getPf().getUnmodifiedForwardPrimer(insert);
        System.out.println("expected: " + expect);
        System.out.println("UnmodifiedFP insertion correct: " + expect.equals(p.seq));
        System.out.println("actual  : " + p.seq);

        Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
        p = pcr.getPf().getUnmodifiedForwardPrimer(mismatch);
        System.out.println("UnmodifiedFP mismatch correct: " + expect.equals(p.seq));
        System.out.println("actual  : " + p.seq);

        Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "del1");
        p = pcr.getPf().getUnmodifiedForwardPrimer(delete);
        System.out.println("UnmodifiedFP delete correct: " + expect.equals(p.seq));
        System.out.println("actual  : " + p.seq);

    }

    public static void testGetModifiedForwardPrimerFromOligo() throws Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        PCR pcr = new PCR(genome);
        Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2, true, "in1");
        String iexpect = "GTTACCTGCCCCCCCAAAAA"; //base 90-99, then the insert
        Primer iprimer = pcr.getPf().getModifiedForwardPrimer(insert);
        System.out.println("Forward from insert oligo: " + iprimer.seq);
        System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.seq.equals(iexpect)));

        Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
        String mexpect = "GTTACCTGCCCCCCCAAAAA"; //80-89, then the new seq
        Primer mprimer = pcr.getPf().getModifiedForwardPrimer(mismatch);
        System.out.println("Forward from mismatch oligo: " + mprimer.seq);
        System.out.println("Test getForwardPrimerFromOligo mismatch : " + (mprimer.seq.equals(mexpect)));

        genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTGGGGGGGGGGAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        pcr = new PCR(genome);
        Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "del1");
        String dexpect = "TTTTTTTTTTTTTTTAAAAA"; //bases 80-89,100-109
        Primer dprimer = pcr.getPf().getModifiedForwardPrimer(delete);
        System.out.println("Forward from deletion oligo: " + dprimer.seq);
        System.out.println("Test getForwardPrimerFromOligo deletion : " + (dprimer.seq.equals(dexpect)));

    }

    //replicate pcr primers used in http://www.sciencemag.org/content/suppl/2011/07/13/333.6040.348.DC1/Isaacs.SOM.pdf
    public static void testPCRbyReplicatingResults() throws Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        PCR pcr = new PCR(genome);
        Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730611, 4730612, 2, true, "TAG->TAA");
        String expectUnmodified = "AACTGGTTGTTAAGCAGTAG";
        String expectModified = "AACTGGTTGTTAAGCAGTAA";
        Primer unmodified = pcr.getPf().getUnmodifiedForwardPrimer(mismatch);
        Primer modified = pcr.getPf().getModifiedForwardPrimer(mismatch);
        System.out.println("Successfully replicated unmodified primer: " + unmodified.seq.equals(expectUnmodified));
        System.out.println("Successfully replicated modified primer: " + modified.seq.equals(expectModified));
    }

    public static void testGeneratePrimerSets() throws IOException, Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        PCR pcr = new PCR(genome);
        ArrayList<Oligo> pool = new ArrayList<Oligo>();
        Oligo insert = Oligo.InsertionFactory(genome, "TTTTTAAAAA", 100, 2, true, "insertion");
        Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730611, 4730612, 2, true, "TAG->TAA");
        Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "deletion");
        pool.add(insert);
        pool.add(mismatch);
        pool.add(delete);
        pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2, true, "oligo1"));
        pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2, true, "oligo2"));
        pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2, true, "oligo3"));
        pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2, true, "oligo4"));
        pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2, true, "oligo5"));
        pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2, true, "oligo6"));
        pool.add(Oligo.InsertionFactory(genome, "GACGATATGAAC", 7900, 2, true, "oligo7"));
        pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2, true, "oligo8"));
        pool.add(Oligo.InsertionFactory(genome, "ATTTGACCT", 15426, 2, true, "oligo9"));
        pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2, true, "oligo10"));

        ArrayList<ArrayList<Primer>> res = pcr.generateAllPrimers(pool);
        //System.out.println(res);

        for (int i = 0; i < res.size(); i++) {
            ArrayList<Primer> arr = res.get(i);
            String s = "";
            for (Primer p : arr) {
                //System.out.println(p);
                if (p != null) {
                    s = s.concat(p.oligo.name + "." + p.amplicon + " ");
                } else {
                    s = s.concat("null ");
                }
            }
            System.out.println(s);
        }
    }

}
