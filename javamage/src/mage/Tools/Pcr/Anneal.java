package mage.Tools.Pcr;

import java.io.IOException;
import java.util.Random;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;

/**
 * Implements a Simulated Annealing algorithm to search the space of all valid
 * PCR primers for the optimum.
 *
 * Simulated annealing utilizes a random walk that gradually reduced the
 * "temperature" parameter to reduce both step size and the probability of
 * accepting a suboptimal step
 *
 * The algorithm: Consider a system with a state described by an N-dimensional
 * vector x , for which the function to be minimized is f(x). The function
 * should return a scalar quantity It may also well be a good idea to do make
 * the transition size ∆x depend on the “temperature”, so that large changes are
 * done at high temperatures, and small ones at lower T’s. 0a Select the initial
 * configuration x. 0b Set the number of Monte Carlo steps nMCS = 0 0c Set the
 * initial temperature to some high value T0 1 Choose a transition ∆x at random
 * 2 Calculate the function value before the transition fb = f(x) 3 Do the trial
 * transtion as x = x + ∆x 4 Calculate the function value after the transition
 * fa = f(x) 5 Calculate ∆f = fa − fb, then : 6 If ∆f ≤ 0 accept the state 7 If
 * ∆f > 0 : 7a Generate a random number u between 0 and 1 7b Accept the state
 * only if e−∆f/T > u 8 If the state is rejected, return to the previous state:
 * x = x − ∆x 9 Reduce the temperature by some small value: T = T − εT 10 Set
 * nMCS = nMCS + 1 12 If T > 0 return to step 1 from Basics of Monte Carlo
 * simulations, Kai Nordlund 2006
 *
 * @author Michael Quintin
 */
/*

 */
public class Anneal {

    //algorithm parameters

    protected float initialT; //starting temp for the algorithm
    protected float finalT = (float) 0.1; //final temp for the algorithm, used by some cooling methods
                                //must be greater than 0. When the difference is one degree
                                //probability of accepting at this temp is e^-10 = 4.5*10^-5
    protected int mcs = 0; //number of Monte Carlo steps completed
    protected int mcsMax = 10000;
    protected int epochLength = 100;
    protected int epoch = 0;
    protected int epochMax = 10;

    public static enum CoolingType {LINEAR, SAWTOOTH, SLOWDOWN};
    protected static CoolingType coolingType = CoolingType.LINEAR; //which method is used to lower t?

    protected Primer[][] matrix; //holds the Primers that have been created while running the search

    protected static Random random = new Random();

    //general primer prameters
    protected double targetTemp; //ideal melting temperature of the optimized primer
    protected String genome; //TODO: consolidate all uses of Genome into one place
    protected int startPos; //5'-most base that can be used
    protected int endPos; //3'-most base that can be used
    protected int lenMin; //minimum acceptible length of a primer
    protected int lenMax; //maximum acceptible length of a primer
    protected double shiftRange; //percent of the amplicon size a reverse
    //primer's start position is allowed to vary by
    protected int maxShift;

    protected PrimerFactory pf;

    //individual primer parameters
    protected Oligo oligo;
    protected int amplicon;
    protected boolean forward;
    protected boolean sense;
    protected boolean modified;

    public Anneal(String genome, double targetTemp, Oligo oligo, int amplicon, double shiftRange, 
            int lenMin, int lenMax, boolean forward, boolean sense, boolean modified) 
            throws IOException {
        this.targetTemp = targetTemp;
        if ((genome == null | "".equals(genome)) & this.genome == null) {
            this.genome = FASTA.readFFN(Oligo.Directory, Oligo.Genome);
        }
        this.pf = new PrimerFactory(genome);
        this.oligo = oligo;
        this.amplicon = amplicon;
        this.forward = forward;
        this.sense = sense;
        this.modified = modified;
        this.lenMin = lenMin;
        this.lenMax = lenMax;
        this.shiftRange = shiftRange;
        
        initializeParams();
    }

    //set up for a run
    protected void initializeParams(){
        maxShift = 0; //how many bases is the start allowed to move?
        if (!forward) {//only reverse primers shift
            Double shiftD = Math.floor(amplicon * shiftRange);
            maxShift = shiftD.intValue();
        }

        //create a matrix to hold Primers so we don't have to recalculate any
        int nRows = 1 + (2 * maxShift);
        int nCols = 1 + (lenMax - lenMin);
        matrix = new Primer[nRows][nCols];

        //set the initial temperature
        //Formally, an "initial acceptance probability" P informs the initial temperature
        //such that P1 is the fraction of uphill transitions that are accepted. A 
        //number of uphill transitions are made and the average increase of the objective
        //function score, Δ, is calculated. The initial temperature is calculated with
        //the equation T = -Δ/(ln P)
        //Here we draw primers at random and use their absolute score difference
        //from each neighbor to estimate delta, and set the initial temperature
        double pAccept = 0.5 ;//arbitrary initial acceptance probability
        double score = 0.0;
        double nToCheck = 5.0;
        double nPairs = 0.0;
        for (int i = 0; i < nToCheck; i++){
            //don't pick Primers on the edge
            int x = (int) (nextInt(maxShift-1) * Math.pow(-1.0, nextInt(2)));//50% to flip to negative 
            int y= nextInt(lenMax-lenMin-2) + lenMin + 1; 
            Primer p1 = getPrimer(x,y);
            int[] adjacency = {-1,1};
            for (int j : adjacency){
                for (int k : adjacency){
                    //one dimensional- y only
                    if ((x == 0 | x == maxShift) & (y != lenMin & y != lenMin)){
                        Primer p2 = getPrimer(x,y+k);
                        score -= Math.abs(score(p1)-score(p2))/2;
                        nPairs += 0.5;
                    }
                    //one dimensional- x only
                    else if ((x != 0 | x != maxShift) & (y == lenMin & y == lenMin)){
                        Primer p2 = getPrimer(x+j,y);
                        score -= Math.abs(score(p1)-score(p2))/2;
                        nPairs += 0.5;
                    }
                    else{
                        Primer p2 = getPrimer(x+j,y+k);
                        score -= Math.abs(score(p1)-score(p2));
                        nPairs++;
                    }
                }
            }
        }
        score = score/nPairs;
        initialT = (float) (score/(Math.log(pAccept)));
    }

    /**
     * Run the algorithm. See the constructor for the Primer class for more
     * details on the parameters
     *
     * @param oligo Oligo this primer corresponds to
     * @param amplicon the size of the amplicon created by this primer, if it's
     * reverse. null if this is a forward primer
     * @param forward is this a forward primer?
     * @param sense which direction is this primer set? if true, the pair
     * follows the given sequence orientation (5' -> 3'). If false, it's flipped
     * (3' -> 5')
     * @param modified is this the modified forward primer?
     * @return
     */
    public Primer getOptimizedPrimer() {
        //set up parameters
        initializeParams();

        float t = initialT;
        mcs = 0;
        int meanLen = (int) (lenMin + lenMax) / 2;
        int[] idx = {0, meanLen}; //index passed to getPrimer() to lookup in the matrix
        Primer current = getPrimer(idx);
        Primer best = getPrimer(idx);

        while ((t > 0) && (epoch < epochMax) && (mcs < mcsMax) && score(current) >= 0.5) {
            //find the next candidate
            int[] move = getMove(t);
            int[] challengerIdx = {idx[0] + move[0], idx[1] + move[1]};

            //bounce off the boundaries
            if (challengerIdx[0] < -maxShift) {
                challengerIdx[0] = -maxShift + (-maxShift - challengerIdx[0]);
            }
            if (challengerIdx[1] < lenMin) {
                challengerIdx[1] = lenMin + (lenMin - challengerIdx[1]);
            }
            if (challengerIdx[0] > maxShift) {
                challengerIdx[0] = maxShift - (challengerIdx[0] - maxShift);
            }
            if (challengerIdx[1] > lenMax) {
                challengerIdx[1] = lenMax - (challengerIdx[1] - lenMax);
            }

            Primer challenger = getPrimer(challengerIdx);

            if (accept(challenger, current, t)) {
                idx = challengerIdx;
                current = challenger;
                if(score(challenger) < score(best)) best = getPrimer(challengerIdx);
            }
            mcs++;
            if (mcs % epochLength == 0) {
                t = cool(t);
            }
        }

        //now run a greedy search, moving only one step at a time
        t = finalT;
        int neighborsChecked = 0;
        while (neighborsChecked < 4) {
            int[] challengerIdx = new int[2];
            challengerIdx[0] = idx[0];
            challengerIdx[1] = idx[1];
            switch (neighborsChecked) {
                case 0:
                    if (idx[0] == -maxShift) {//already at the edge
                        break;
                    }
                    challengerIdx[0] = challengerIdx[0] - 1;
                    break;
                case 1:
                    if (idx[0] == -maxShift) {
                        break;
                    }
                    challengerIdx[0] = challengerIdx[0] + 1;
                    break;
                case 2:
                    if (idx[1] == lenMin) {
                        break;
                    }
                    challengerIdx[1] = challengerIdx[1] - 1;
                    break;
                case 3:
                    if (idx[1] == lenMax) {
                        break;
                    }
                    challengerIdx[1] = challengerIdx[1] + 1;
                    break;
            }
            neighborsChecked += 1;
            Primer challenger = getPrimer(challengerIdx);
            if (idx[0] != challengerIdx[0] || idx[1] != challengerIdx[1]){
                if (accept(challenger, current, t)){
                    current = challenger;
                    idx = challengerIdx;
                    neighborsChecked = 0;
                    if(score(challenger) < score(best)) {
                        best = getPrimer(challengerIdx);
                    }
                }
            }
        }
        /*System.out.println("Merlin-DEBUG: Generated Primer for " + oligo.name + 
                " Amplicon " + amplicon + " Sense=" + best.sense + 
                " Modified=" + best.modified + " with MT " + best.getMt() + 
                ": " + best.seq);*/
        return best;
    }

    /**get the relative coordinates of the next primer to check through a random walk
     * where the number of steps is controlled by the temperature. The initial distance
     * is just enough to go from the center to a corner, the distance in the final
     * epoch is a single step
     * 
     * @param t
     * @return 
     */
    protected int[] getMove(float t){
        int[] move = {0,0};
        int initDist = matrix.length + matrix[0].length -1;
        if (initDist == 1) return move; //only a single primer is possible
        
        //scales from initDist+1 at t=initialT to 1 when t=finalT
        int nSteps = (int) Math.round(initDist - ((initDist-1) * (initialT - t)/(initialT-finalT)));
        while (nSteps > 0){
            int r = nextInt(4);
            if (matrix[0].length == 1) r = nextInt(2);
            if (matrix.length == 1) r = nextInt(2) + 2;
            switch (r){
                case 0:
                    move[0] = move[0] + 1;
                    break;
                case 1:
                    move[0] = move[0] - 1;
                    break;
                case 2:
                    move[1] = move[1] + 1;
                    break;
                case 3:
                    move[1] = move[1] - 1;
                    break;
            }
            nSteps--;
        }
        return move;
    }

    /**
     * Keep the first Primer ("challenger") over the second ("incumbent")? If ∆f
     * is less than 0, keep the challenger Otherwise make a random number
     * between 0 and 1 and keep the challenger if it's smaller than e^(-∆f/t)
     *
     * @param challenger
     * @param incumbent
     * @param t current algorithm temperature
     * @return
     */
    public boolean accept(Primer challenger, Primer incumbent, float t) {
        Double diff = score(challenger) - score(incumbent);
        if (diff < 0) {
            return true;
        } 
        else {
            double rand = Math.random();
            double odds = Math.exp((-1 * diff) / t);
            //a temp difference of 0.03 degrees replaces about at about 75%
            //smaller temp differences may make the greedy optimization run too long
            odds = Math.min(odds, Math.exp((-1 * 0.03) / t));  
            return rand < odds;
        }
    }

    /**
     * Get a Primer from the storage matrix according to its shift and length,
     * instead of having to track indices. If the primer doesn't exist, create
     * and store it
     *
     * @param shift
     * @param len
     * @return
     */
    protected Primer getPrimer(int shift, int len) {
        int sidx = shift + maxShift;
        int lenidx = len - lenMin;
        Primer p = null;
        if (sidx < matrix.length){ //watch out for null pointer exception
            if (lenidx < matrix[sidx].length){
                p = matrix[sidx][lenidx];
            }
        }
        if (p != null) {
            return p;
        } else {
            p = pf.getPrimer(oligo, amplicon, forward, sense, modified, shift, len);
            matrix[sidx][lenidx] = p;
        }
        return p;
    }

    //helper to not have to unpack an index array
    protected Primer getPrimer(int[] arr) {
        return getPrimer(arr[0], arr[1]);
    }

    /**
     * Score the given primer according to the given target temperature A LOW
     * score is better- it's based on the absolute distance from the target
     *
     * @param p
     * @return
     */
    public double score(Primer p) {
        double primerTemp = p.getMt();
        //return Math.pow(targetTemp - primerTemp,2);
        return Math.abs(targetTemp - primerTemp);
    }

    /**
     * reduce the annealing temperature
     *
     */
    protected float cool(float t) {
        epoch += 1;
        float newT = t;
        switch (coolingType) {
            case LINEAR:
                float delta = (initialT - finalT)/epochMax;
                newT -= delta;
                break;
            case SAWTOOTH://TODO

                break;
            case SLOWDOWN:
                //A cooling function which gives slower cooling at higher epochs: 
                //T(k)=T(k-1)/(1+bT(k-1)) where b>0. Set b=[T(0)-T(M)]/[(M-1)*T(0)*T(M)]
                float b = (initialT - finalT) / ((epochMax - 1) * initialT * finalT);
                newT = (t / (1 + (b * t)));
                break;
        }
        return newT;
    }
    
    //Random.nextInt requires the given integer to be positive. This function
    //helps in the case where the primer matrix is one dimensional
    protected int nextInt(int i){
        if (0 >= i) return 0;
        else return Anneal.random.nextInt(i);
    }

    //-------Getters and Setters
    public float getInitialT() {
        return initialT;
    }

    public void setInitialT(float initialT) {
        this.initialT = initialT;
    }

    public int getMcs() {
        return mcs;
    }

    public void setMcs(int mcs) {
        this.mcs = mcs;
    }

    public int getMcsMax() {
        return mcsMax;
    }

    public void setMcsMax(int mcsMax) {
        this.mcsMax = mcsMax;
    }

    public Primer[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(Primer[][] matrix) {
        this.matrix = matrix;
    }

    public double getTargetTemp() {
        return targetTemp;
    }

    public void setTargetTemp(double targetTemp) {
        this.targetTemp = targetTemp;
    }

    public String getGenome() {
        return genome;
    }

    public void setGenome(String genome) {
        this.genome = genome;
    }

    public int getStartPos() {
        return startPos;
    }

    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    public int getEndPos() {
        return endPos;
    }

    public void setEndPos(int endPos) {
        this.endPos = endPos;
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

    public double getShiftRange() {
        return shiftRange;
    }

    public void setShiftRange(double shiftRange) {
        this.shiftRange = shiftRange;
    }

    public int getMaxShift() {
        return maxShift;
    }

    public void setMaxShift(int maxShift) {
        this.maxShift = maxShift;
    }

    public Oligo getOligo() {
        return oligo;
    }

    public void setOligo(Oligo oligo) {
        this.oligo = oligo;
    }

    public int getAmplicon() {
        return amplicon;
    }

    public void setAmplicon(int amplicon) {
        this.amplicon = amplicon;
    }

    public boolean isForward() {
        return forward;
    }

    public void setForward(boolean forward) {
        this.forward = forward;
    }

    public boolean isSense() {
        return sense;
    }

    public void setSense(boolean sense) {
        this.sense = sense;
    }

    public boolean isModified() {
        return modified;
    }

    public void setModified(boolean modified) {
        this.modified = modified;
    }

    public void setEpoch(int epoch) {
        this.epoch = epoch;
    }

    public int getEpoch() {
        return this.epoch;
    }

    public void setEpochMax(int epochMax) {
        this.epochMax = epochMax;
    }

    public int getEpochMax() {
        return this.epochMax;
    }
    
    public float getFinalT() {
        return finalT;
    }

    public void setFinalT(float finalT) {
        this.finalT = finalT;
    }

    public static Random getRandom() {
        return random;
    }

    public static void setRandom(Random random) {
        Anneal.random = random;
    }

    public PrimerFactory getPf() {
        return pf;
    }

    public void setPf(PrimerFactory pf) {
        this.pf = pf;
    }
    
    public static CoolingType getCoolingType() {
        return coolingType;
    }

    public static void setCoolingType(CoolingType coolingType) {
        Anneal.coolingType = coolingType;
    }

}
