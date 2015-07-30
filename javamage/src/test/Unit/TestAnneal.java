/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test.Unit;

import java.io.IOException;
import java.util.ArrayList;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;
import mage.Tools.Pcr.Anneal;
import test.Constants;
/**
 *
 * @author mquintin
 */
public class TestAnneal extends Anneal {
    public Integer[][] epochTracker; //recored the epoch when a primer gets made
    public Integer[][] mcsTracker;
    public Double[][] scoreTracker;

    public static void main(String[] args) throws Exception {
        String gen = mage.Tools.FASTA.readFFN(Constants.blastdirectory, "genome.ffn");
        Oligo o = Oligo.InsertionFactory(gen, "TAG", 3090358, 2, true, "endA");
        //Oligo o = Oligo.MismatchFactory(genome, "A", left_position, right_position, replichore, true, genome)

        //test Move
        if (false) {
            TestAnneal ta = new TestAnneal(60.0, o, 400, .05, 16, 30, false, true, false);
            System.out.println("TestMove in Two dimensions, T = " + ta.initialT);
            ta.testGetMove(ta.initialT);

            float cooledT = ta.cool(ta.initialT);
            cooledT = ta.cool(cooledT);
            System.out.println("Partially cooled, T = " + cooledT);
            ta.testGetMove(cooledT);
            System.out.println("Fully Cooled, T = " + ta.finalT);
            ta.testGetMove(ta.finalT);
            ta = new TestAnneal(60.0, o, 400, 0, 16, 30, false, true, false);
            System.out.println("TestMove, No shifting");
            ta.testGetMove(ta.initialT); //no shifting allowed
            ta = new TestAnneal(60.0, o, 400, .05, 20, 20, false, true, false);
            System.out.println("TestMove, Fixed length");
            ta.testGetMove(ta.initialT); //fixed length
        }

        //test Cool
        if (false) {
            //Slowdown
            Anneal.setCoolingType(CoolingType.SLOWDOWN);
            TestAnneal ta = new TestAnneal(60.0, o, 400, .05, 16, 30, false, true, false);
            ta.testCool();
            //Linear
            Anneal.setCoolingType(CoolingType.LINEAR);
            ta = new TestAnneal(60.0, o, 400, .05, 16, 30, false, true, false);
            ta.testCool();

        }

        //test getOptimizedPrimer
        if (false) {
            TestAnneal ta = new TestAnneal(60.0, o, 400, .05, 10, 30, false, true, false);
            Anneal.random.setSeed(9999l);
            ta.testOptimize();
            ta = new TestAnneal(60.0, o, 400, .05, 10, 30, false, true, false);
            Anneal.random.setSeed(8888l);
            ta.testOptimize();
            ta = new TestAnneal(60.0, o, 400, .05, 10, 30, false, true, false);
            Anneal.random.setSeed(7777l);
            ta.testOptimize();
        }
        
        //test accept
        if (false){
            TestAnneal ta = new TestAnneal(60.0, o, 400, .05, 10, 30, false, true, false);
            System.out.println("Test Accept at initial temp");
            Anneal.random.setSeed(999l);
            ta.testAccept(ta.initialT);
            ta.testAccept2(ta.initialT);
            System.out.println("Test Accept at final temp");
            Anneal.random.setSeed(888l);
            //ta.testAccept(ta.finalT);
            
        }
        
        //test initial temp
        if (false){
            TestAnneal ta = new TestAnneal(60.0, o, 400, .05, 10, 30, false, true, false);
            ta.testInitialTemp();
            ta.testInitialTemp();
            ta.testInitialTemp();
        }
        
        //test a particular problem primer
        if (true) {
            String directory = Constants.blastdirectory + "nat/";
            String genome = FASTA.readFFN(directory, "genome.FASTA");
            Oligo.Directory = directory;
            Oligo.Genome = "genome.FASTA";
            o = Oligo.MismatchFactory(genome, "A", 4572079, 4572079, 1, true, "avrII07");
            TestAnneal ta = new TestAnneal(60.0, o, 500, .05, 10, 30, false, true, true);
            ta.testOptimize();
        }

    }

    public TestAnneal(double targetTemp, Oligo oligo, int amplicon, double shiftRange, int lenMin, int lenMax, boolean forward, boolean sense, boolean modified) throws IOException {
        super(mage.Tools.FASTA.readFFN(Constants.blastdirectory, "genome.ffn"), targetTemp, oligo, amplicon, shiftRange, lenMin, lenMax, forward, sense, modified);
        //creationTracker = new int[matrix.length][matrix[0].length];
        epochTracker = new Integer[matrix.length][matrix[0].length];
        mcsTracker = new Integer[matrix.length][matrix[0].length];
        scoreTracker = new Double[matrix.length][matrix[0].length];
    }

    public void testGetMove(float t) throws IOException {
        System.out.println("Matrix dims: " + matrix.length + " X " + matrix[0].length);
        double deltaSum = 0.0;
        for (int i = 0; i < 30; i++) {
            int[] arr = getMove(t);
            String s = "TestMove   " + i + ":\t";
            s = s + arr[0] + "\t" + arr[1];
            double delta = (Math.abs(arr[0]) + Math.abs(arr[1]));
            deltaSum += delta;
            s = s + "\tDelta = " + delta;
            System.out.println(s);
        }
        System.out.println("Average Delta = " + (deltaSum / 30.0));
    }

    public void testCool() {
        System.out.println("Test Cooling Type: " + getCoolingType().toString());
        int n = 10;
        ArrayList<Float> temps = new ArrayList();
        temps.add(initialT);
        ArrayList<Float> deltas = new ArrayList();
        float latest = initialT;
        for (int i = 1; i < n; i++) {
            latest = cool(latest);
            temps.add(latest);
            deltas.add(temps.get(i) - temps.get(i - 1));
        }
        System.out.println("Temps: " + temps);
        System.out.println("Deltas: " + deltas);
    }

    public void testOptimize() {
        long starttime = System.nanoTime();
        Primer opt = getOptimizedPrimer();
        long endtime = System.nanoTime();

        //compare to all other calculated primers
        double best = 999.0;
        Primer bestPrimer = null;
        int n = 0;
        for (Primer[] row : matrix) {
            for (Primer p : row) {
                if (p != null) {
                    n++;
                    double s = score(p);
                    best = Math.min(s, best);
                    bestPrimer = p;
                }
            }
        }
        System.out.println("Optimize finished in t = " + (endtime - starttime));
        System.out.println(mcs + " Monte Carlo steps were performed");
        System.out.println(n + " primers were created of a possible " + 
                (matrix.length * matrix[0].length));
        System.out.println("Found primer's score was " + score(opt) + ". Best seen was " + best);
        System.out.println("Found primer's coordinates are Shift: " +
                (opt.genome_start - opt.oligo.target_position - opt.amplicon) + ", Len: " +
                opt.seq.length());
        System.out.println("Optimal primer returned: " + (score(opt) == best));
        System.out.println("Returned primer seq: " + opt.seq);
        System.out.println("Best primer seq: " + bestPrimer.seq);
        System.out.println();
        /*System.out.println("Matrix population: ");
        for (Primer[] row : matrix){
            String s = "";
            for (Primer p : row){
                int i = 0;
                if (p != null) i = 1;
                s = s + i + "\t";
            }
            System.out.println(s);
        }*/
        System.out.println("Matrix population history: ");
        ArrayList<Number[][]> trackers = new ArrayList();
        trackers.add(mcsTracker);
        trackers.add(epochTracker);
        trackers.add(scoreTracker);
        for (Number[][] tracker  : trackers) {
            for (int r = 0; r < tracker.length; r++) {
                Number[] row = tracker[r];
                String s = r - maxShift + ":\t";
                for (Number i : row) {
                    if (i == null) {
                        s = s + "\t";
                    } else {
                        s = s + i + "\t";
                    }
                }
                System.out.println(s);
            }
            String yaxis = "\t";
            for (int i = lenMin; i <= lenMax; i++) {
                yaxis = yaxis + i + "\t";
            }
            System.out.println(yaxis);
        }
    }
    
    public void testAccept(float temp){
        double trials = 100.0;
        Primer[] primers = new Primer[3];
        primers[0] = getPrimer(0,20);
        primers[1] = getPrimer(5,25);
        primers[2] = getPrimer(-5,15);
        for (int i = 0; i < primers.length; i++){
            for (int j = 0; j < primers.length; j++){
                if (i != j){
                    double diff = score(primers[i]) - score(primers[j]);
                    double accepted = 0.0;
                    for (int t = 0; t < trials; t++){
                        if (accept(primers[i],primers[j], temp)) accepted++;
                    }
                    System.out.println("P1: " + score(primers[i]) + "\tP2: " + 
                            score(primers[j]) + "\t %Accepted: " + (accepted/trials));
                }
            }
        }
    }
    
    public void testAccept2(float temp){
        double trials = 100.0;
        double accepted = 0.0;
        for (int i = 0; i < trials; i++){
            Primer p1 = getPrimer(nextInt(maxShift),random.nextInt(lenMax-lenMin) + lenMin);
            Primer p2 = getPrimer(nextInt(maxShift),random.nextInt(lenMax-lenMin) + lenMin);
            if (score(p1)>score(p2)){
                if (accept(p1,p2,temp)) accepted++;
            }
            else if (accept(p2,p1,temp)) accepted++;
        }
        System.out.println("% Accepted = " + (accepted/trials));
    }
    
    public void testInitialTemp(){
        ArrayList<Float> arr = new ArrayList();
        int trials = 10;
        float total = 0f;
        for (int i = 0; i < trials ; i++){
            initializeParams();
            arr.add(initialT);
            total += initialT;
        }
        System.out.println("Initial Temps: " + arr.toString());
        System.out.println("Average Initial Temp: " + (total/trials));
    }
    
    @Override
    protected Primer getPrimer(int[] arr){
        int sidx = arr[0] + maxShift;
        int lenidx = arr[1] - lenMin;
        epochTracker[sidx][lenidx] = epoch;
        mcsTracker[sidx][lenidx] = mcs;
        Primer p = getPrimer(arr[0],arr[1]);
        scoreTracker[sidx][lenidx] = score(p);
        return p;
    }
   
}
