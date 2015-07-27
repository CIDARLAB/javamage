package test;

import java.io.IOException;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;
import mage.Tools.Pcr.Anneal;

/** Class to examine performance of the annealing algorithm to determine ideal
 * parameters.
 *
 * For further reading, see Park & Kim 1998, Computers Ops Res. Vol 25
 * 
 * @author mquintin
 */
public class AnnealParameters {
    public static String genome;
    public static float t_init; //initial temp
    public static float t_final; //final temp
    public static int n_epochs; //how many temp steps do we use?
    public static int len_epoch; //how many trials in each epoch
    
    public static void main(String[] args) throws Exception{
        genome = FASTA.readFFN(Constants.blastdirectory, "genome.ffn");
        Oligo o = Oligo.InsertionFactory(genome, "TAG", 3090358, 2, true, "endA");
        proportionalCooling(o);
        //idealTemp(o);
    }
    
    /**An "initial acceptance probability" P informs the initial temperature
     * such that P1 is the fraction of uphill transitions that are accepted. A 
     * number of uphill transitions are made and the average increase of the objective
     * function score, Δ, is calculated. The initial temperature is calculated with
     * the equation T = Δ/(ln P)
     * 
     * @throws Exception 
     */
    public static void idealTemp(Oligo o){
                
        
    }
    
    /**In the proportional cooling scheme, T(k) = a*T(k-1) where a is in the range (0,1).
     * A method to create a depends on the initial temp T(0), final temp T(M) 
     * and the number of epochs M such that a = [T(M)/T(0)]^[1/(M-1)] 
     * 
     * @param o 
     */
    public static void proportionalCooling(Oligo o){
        
    }
    
    /**A cooling function which gives slower cooling at higher epochs: 
     * T(k)=T(k-1)/(1+bT(k-1)) where b>0. Set b=[T(0)-T(M)]/[(M-1)*T(0)*T(M)]  
     * 
     * @param o 
     */
    public static void deceleratingCooling(Oligo o) throws IOException{
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        Anneal anneal = new Anneal(genome, 60.0, o, 400, 0.05, 16, 30, false, true, false);
        anneal.setLenMin(16);
        anneal.setLenMax(30);
        anneal.setShiftRange(0.1);
        anneal.setMaxShift(20);
        anneal.setEpochMax(1);
        float b;
        Primer p = null;
        for (int i = 0; i < n_epochs; i++){
            p = anneal.getOptimizedPrimer();
            b = (t_init - t_final)/((n_epochs - 1) * t_init * t_final);
            float t_current = anneal.getInitialT();
            anneal.setInitialT(t_current/(1 + (b * t_current)));
        }
        if (p != null) System.out.println(p.seq + ", mt: " + p.getMt());
    }
}
