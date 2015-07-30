/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.io.IOException;
import java.util.ArrayList;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;
import mage.Tools.OutputTools;
import mage.Tools.Pcr.Melt;
import mage.Tools.Pcr.PCR;
import mage.Tools.Pcr.PrimerFactory;

/**
 *
 * @author mquintin
 */
public class TestingSandbox {

    public static String genome;
    public static PCR pcr;

    public static void main(String[] args) throws IOException, Exception {
        genome = FASTA.readFFN(Constants.blastdirectory, "genome.ffn");
        //pcr = new PCR(genome);

        //runtimeToCalculateAllPrimers();
        //System.out.println("createAvrTry7PrimersInclusive");
        //createAvrTry7PrimersInclusive();
        System.out.println("createAvrTry7PrimersExclusive");
        createAvrTry7PrimersExclusive();
        /*System.out.println("createAvrTry6PrimersInclusive");
        createAvrTry6PrimersInclusive();
        System.out.println("createAvrTry6PrimersExclusive");
        createAvrTry6PrimersExclusive();
        */
        //testAvrII06Oligo();
    }
    
    //this oligo doesn't seem to differentiate from the genome
    public static void testAvrII06Oligo() throws Exception{
        String directory = Constants.blastdirectory + "nat/";
        String genome = FASTA.readFFN(directory, "genome.FASTA");
        Oligo.Directory = directory;
        Oligo.Genome = "genome.FASTA";
        Oligo o = Oligo.MismatchFactory(genome, "N", 3795826, 3795827, 2, true, "avrII06");
        ArrayList<Primer> results = new ArrayList();
        PrimerFactory pf = new PrimerFactory(genome);
        results.add(pf.getUnmodifiedForwardPrimer(o));
        results.add(pf.getModifiedForwardPrimer(o));
        results.add(pf.getUnmodifiedAntisenseForwardPrimer(o));
        results.add(pf.getModifiedAntisenseForwardPrimer(o));

        for (Primer p : results) {
            System.out.println(p.oligo.name + "\tSense:" + p.sense
                    + "\tModified:" + p.modified + "\tReversed:"
                    + p.oligo.sequence.substring(0, 1).equals(p.seq.substring(0, 1))
                    + "\tCorrect:" + !p.seq.substring(p.seq.length() - 1).equals("N")
                    + "\t" + p.seq);
        }
    }

    public static void createAvrTry7PrimersInclusive() throws Exception {
        String directory = Constants.blastdirectory + "nat/";
        genome = FASTA.readFFN(directory, "genome.FASTA");
        Oligo.Directory = directory;
        Oligo.Genome = "genome.FASTA";
        ArrayList<Oligo> pool = new ArrayList();
        pool.add(Oligo.MismatchFactory(genome, "G", 168926, 168926, 1, true, "avrII01"));
        pool.add(Oligo.MismatchFactory(genome, "A", 292079, 292079, 1, true, "avrII02"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1196072, 1196072, 1, true, "avrII03"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1432189, 1432189, 1, true, "avrII04"));
        pool.add(Oligo.MismatchFactory(genome, "T", 1631155, 1631155, 2, true, "avrII05"));
        pool.add(Oligo.MismatchFactory(genome, "C", 3795826, 3795826, 2, true, "avrII06"));
        pool.add(Oligo.MismatchFactory(genome, "A", 4572079, 4572079, 1, true, "avrII07"));

        OutputTools.generateMASCPCRFile(pool, directory + "AvrTry7Inclusive.txt");
    }

    public static void createAvrTry7PrimersExclusive() throws Exception {
        String directory = Constants.blastdirectory + "nat/";
        genome = FASTA.readFFN(directory, "genome.FASTA");
        Oligo.Directory = directory;
        Oligo.Genome = "genome.FASTA";
        ArrayList<Oligo> pool = new ArrayList();
        pool.add(Oligo.MismatchFactory(genome, "G", 168926, 168927, 1, true, "avrII01"));
        pool.add(Oligo.MismatchFactory(genome, "A", 292079, 292080, 1, true, "avrII02"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1196072, 1196073, 1, true, "avrII03"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1432189, 1432190, 1, true, "avrII04"));
        pool.add(Oligo.MismatchFactory(genome, "T", 1631155, 1631156, 2, true, "avrII05"));
        pool.add(Oligo.MismatchFactory(genome, "C", 3795826, 3795827, 2, true, "avrII06"));
        pool.add(Oligo.MismatchFactory(genome, "A", 4572079, 4572080, 1, true, "avrII07"));

        OutputTools.generateMASCPCRFile(pool, directory + "AvrTry7Exclusive.txt");
    }
    
    public static void createAvrTry6PrimersInclusive() throws Exception{
        String directory = Constants.blastdirectory + "nat/";
        genome = FASTA.readFFN(directory, "genome.FASTA");
        Oligo.Directory = directory;
        Oligo.Genome = "genome.FASTA";
        ArrayList<Oligo> pool = new ArrayList();
        pool.add(Oligo.MismatchFactory(genome, "G", 168926, 168927, 1, true, "avrII01"));
        pool.add(Oligo.MismatchFactory(genome, "A", 292079, 292079, 1, true, "avrII02"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1196072, 1196072, 1, true, "avrII03"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1432189, 1432189, 1, true, "avrII04"));
        pool.add(Oligo.MismatchFactory(genome, "T", 1631155, 1631155, 2, true, "avrII05"));
        pool.add(Oligo.MismatchFactory(genome, "C", 3795826, 3795826, 2, true, "avrII06"));
        pool.add(Oligo.MismatchFactory(genome, "A", 4572079, 4572079, 1, true, "avrII07"));
        pool.add(Oligo.MismatchFactory(genome, "TGGTGGGGTAACGGCTCACCA", 224023, 224023, 1, true, "avrIIrrsABCEH"));
        pool.add(Oligo.MismatchFactory(genome, "C", 2727404, 2727404, 2, true, "avrIIgltW"));
        pool.add(Oligo.MismatchFactory(genome, "G", 3941520, 3941520, 1, true, "avrIIgltU"));
        pool.add(Oligo.MismatchFactory(genome, "G", 4166457, 4166457, 1, true, "avrIIgltT"));
        pool.add(Oligo.MismatchFactory(genome, "G", 4207859, 4207859, 1, true, "avrIIgltV"));
        
        OutputTools.generateMASCPCRFile(pool, directory + "AvrTry6Inclusive.txt");        
        
    }
    
        public static void createAvrTry6PrimersExclusive() throws Exception{
        String directory = Constants.blastdirectory + "nat/";
        genome = FASTA.readFFN(directory, "genome.FASTA");
        Oligo.Directory = directory;
        Oligo.Genome = "genome.FASTA";
        ArrayList<Oligo> pool = new ArrayList();
        pool.add(Oligo.MismatchFactory(genome, "G", 168926, 168928, 1, true, "avrII01"));
        pool.add(Oligo.MismatchFactory(genome, "A", 292079, 292080, 1, true, "avrII02"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1196072, 1196073, 1, true, "avrII03"));
        pool.add(Oligo.MismatchFactory(genome, "A", 1432189, 1432190, 1, true, "avrII04"));
        pool.add(Oligo.MismatchFactory(genome, "T", 1631155, 1631156, 2, true, "avrII05"));
        pool.add(Oligo.MismatchFactory(genome, "C", 3795826, 3795827, 2, true, "avrII06"));
        pool.add(Oligo.MismatchFactory(genome, "A", 4572079, 4572080, 1, true, "avrII07"));
        pool.add(Oligo.MismatchFactory(genome, "TGGTGGGGTAACGGCTCACCA", 224023, 224024, 1, true, "avrIIrrsABCEH"));
        pool.add(Oligo.MismatchFactory(genome, "C", 2727404, 2727405, 2, true, "avrIIgltW"));
        pool.add(Oligo.MismatchFactory(genome, "G", 3941520, 3941521, 1, true, "avrIIgltU"));
        pool.add(Oligo.MismatchFactory(genome, "G", 4166457, 4166458, 1, true, "avrIIgltT"));
        pool.add(Oligo.MismatchFactory(genome, "G", 4207859, 4207860, 1, true, "avrIIgltV"));
        
        OutputTools.generateMASCPCRFile(pool, directory + "AvrTry6Exclusive.txt");        
        
    }

    //recreate a primer from a real data set that failed
    public static void recreateLowTempPrimer() throws Exception {
        Oligo o = Oligo.MismatchFactory(genome, "A", 3090358, 3090358, 2, true, "endA");
        Melt.setCmdArg("dntps", "0.8");
        Melt.setCmdArg("na", "50");
        Melt.setCmdArg("mg", "3");

        //recreate a primer from a use case 
        /*Primer fwdA = pcr.getModifiedAntisenseForwardPrimer(o);
         Primer fwdS = pcr.getModifiedForwardPrimer(o);
         Primer up = pcr.getUpstreamReversePrimer(o, o.getGenomeStart(), 200);
         Primer down = pcr.getDownstreamReversePrimer(o, o.getGenomeStart(), 200);
         Melt.setCmdArg("dntps", "0.8");
         Melt.setCmdArg("na", "50");
         Melt.setCmdArg("mg", "3");
         up.getMt();
         down.getMt();
         ArrayList<Primer> primers = new ArrayList();
         //primers.add(up);
         //primers.add(down);
         for (int i = 0; i < 40; i++ ){
         down = pcr.getDownstreamReversePrimer(o, o.getGenomeStart() - 20 + i, 200);
         down = pcr.optimizePrimer(down);
         primers.add(down);
         }
        
         int id = 0;
         for (Primer p : primers){
         id++;
         System.out.println("ID:" + id + " S:" + p.start + " L:" + p.seq.length() + " Mt:" + p.getMt());
         }

         System.out.println(up);*/
    }

    //test how long it takes to calculate all PCR primers at a specific site
    public static void runtimeToCalculateAllPrimers() throws IOException, Exception {
        Oligo o = Oligo.MismatchFactory(genome, "A", 3090358, 3090358, 2, true, "endA");
        ArrayList<Primer> primers = new ArrayList(960);
        //double tempShift = pcr.getPf().getMaxshift();
        /*int tempLen = pcr.getPrimerlength();
         pcr.setMaxshift(100.0);
         long start = System.nanoTime();
         //40 base start window
         //16-40 base length window
         //960 primers
         for (int i = 0; i < 40; i++) {
         for (int j = 16; j <= 40; j++) {
         pcr.setPrimerlength(j);
         Primer p = pcr.getUpstreamReversePrimer(o, o.getGenomeStart() + i, 200);
         primers.add(p);
         }
         }
         //get MTs
         Melt.setMTs(primers);
         long end = System.nanoTime();
         System.out.println("Generating all primers took " + (end-start) + "ns");
        
         //reset PCR
         pcr.setMaxshift(tempShift);
         pcr.setPrimerlength(tempLen);*/
    }
}
