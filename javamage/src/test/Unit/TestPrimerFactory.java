/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test.Unit;

import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.FASTA;
import mage.Tools.Pcr.PCR;
import mage.Tools.Pcr.PrimerFactory;
import mage.Tools.SequenceTools;
import test.Constants;

/**
 *
 * @author mquintin
 */
public class TestPrimerFactory {

    public static void main(String[] args) throws Exception {
        //testGetModifiedForwardPrimerFromOligo();
        //testGetAntisenseModifiedForwardPrimerFromOligo();
        replicateResults();
    }

    public static void testGetModifiedForwardPrimerFromOligo() throws Exception {
        String genome = mage.Tools.FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        //PCR pcr = new PCR(genome);
        PrimerFactory pf = new PrimerFactory(genome);
        Oligo insert = Oligo.InsertionFactory(genome, "TTTTTAAAAA", 90, 2, true, "in1");
        String iexpect = "GTTACCTGCCTTTTTAAAAA"; //base 90-99, then the insert
        Primer iprimer = pf.getModifiedForwardPrimer(insert);
        System.out.println("Forward from insert oligo: " + iprimer.seq);
        int beginIndex = insert.target_position;
        int endIndex = beginIndex + insert.target_length;
        iexpect = insert.sequence.substring(endIndex-pf.getPrimerLength(),endIndex);
        System.out.println("Expect  from insert oligo: " + iexpect);
        System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.seq.equals(iexpect)));

        Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
        String mexpect = "GTTACCTGCCCCCCCAAAAA"; //80-89, then the new seq
        Primer mprimer = pf.getModifiedForwardPrimer(mismatch);
        System.out.println("Forward from mismatch oligo: " + mprimer.seq);
        beginIndex = mismatch.target_position;
        endIndex = beginIndex + mismatch.target_length;
        mexpect = mismatch.sequence.substring(endIndex-pf.getPrimerLength(),endIndex);
        System.out.println("Expect  from mismatch oligo: " + mexpect);
        System.out.println("Test getForwardPrimerFromOligo mismatch : " + (mprimer.seq.equals(mexpect)));

        genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTGGGGGGGGGGAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        Oligo delete = Oligo.DeletionFactory(genome, 90, 99, 2, "del1");
        String dexpect = "TTTTTTTTTTTTTTTAAAAA"; //bases 80-89,100-109
        Primer dprimer = pf.getModifiedForwardPrimer(delete);
        System.out.println("Forward from deletion oligo: " + dprimer.seq);
        System.out.println("Expect  from deletion oligo: " + dexpect);
        System.out.println("Test getForwardPrimerFromOligo deletion : " + (dprimer.seq.equals(dexpect)));

    }
    public static void testGetAntisenseModifiedForwardPrimerFromOligo() throws Exception {
        String genome = mage.Tools.FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        //PCR pcr = new PCR(genome);
        PrimerFactory pf = new PrimerFactory(genome);
        Oligo insert = Oligo.InsertionFactory(genome, "CCCCCAAAAA", 90, 2, true, "in1");
        String iexpect = ""; //base 90-99, then the insert
        Primer iprimer = pf.getModifiedAntisenseForwardPrimer(insert);
        System.out.println("Forward from insert oligo: " + iprimer.seq);
        int beginIndex = insert.target_position + insert.target_length - pf.getPrimerLength();
        int endIndex = beginIndex + pf.getPrimerLength();
        iexpect = insert.sequence.substring(beginIndex, endIndex);
        iexpect = SequenceTools.ReverseCompliment(iexpect);
        System.out.println(insert.sequence);
        System.out.println(SequenceTools.ReverseCompliment(insert.sequence));
        System.out.println("Expect  from insert oligo: " + iexpect);
        System.out.println("Test getForwardPrimerFromOligo insert : " + (iprimer.seq.equals(iexpect)));

        Oligo mismatch = Oligo.MismatchFactory(genome, "CCCCCAAAAA", 90, 100, 2, true, "mis1");
        String mexpect = ""; //80-89, then the new seq
        Primer mprimer = pf.getModifiedAntisenseForwardPrimer(mismatch);
        beginIndex = mismatch.target_position + mismatch.target_length - pf.getPrimerLength();
        mexpect = mismatch.sequence.substring(beginIndex,beginIndex + pf.getPrimerLength());
        mexpect = SequenceTools.ReverseCompliment(mexpect);
        System.out.println(mismatch.sequence);
        System.out.println(SequenceTools.ReverseCompliment(mismatch.sequence));
        System.out.println("Forward from mismatch oligo: " + mprimer.seq);
        System.out.println("Expect  from mismatch oligo: " + mexpect);
        System.out.println("Test getForwardPrimerFromOligo mismatch : " + (mprimer.seq.equals(mexpect)));

        genome = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTGGGGGGGGGGAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        Oligo delete = Oligo.DeletionFactory(genome, 90, 99, 2, "del1");
        String dexpect = "GGGGGGGGGGTTTTTAAAAA"; 
        Primer dprimer = pf.getModifiedAntisenseForwardPrimer(delete);
        System.out.println("Forward from deletion oligo: " + dprimer.seq);
        System.out.println("Expect  from deletion oligo: " + dexpect);
        System.out.println("Test getForwardPrimerFromOligo deletion : " + (dprimer.seq.equals(dexpect)));

    }
    
        //replicate pcr primers used in http://www.sciencemag.org/content/suppl/2011/07/13/333.6040.348.DC1/Isaacs.SOM.pdf
    public static void replicateResults() throws Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        PrimerFactory pf = new PrimerFactory(genome,26);
        Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730612, 4730612, 2, true, "TAG->TAA");
        System.out.println(mismatch.sequence);
        System.out.println(SequenceTools.ReverseCompliment(mismatch.sequence));
        String expectUnmodified = "ATCTGAAACTGGTTGTTAAGCAGTAG";
        String expectModified = "ATCTGAAACTGGTTGTTAAGCAGTAA";
        Primer unmodified = pf.getUnmodifiedForwardPrimer(mismatch);
        Primer modified = pf.getModifiedForwardPrimer(mismatch);
        System.out.println("Expected  Unmodified: " + expectUnmodified);
        System.out.println("Generated Unmodified: " + unmodified.seq);
        System.out.println();
        System.out.println("Expected  Modified: " + expectModified);
        System.out.println("Generated Modified: " + modified.seq);
        System.out.println();
        System.out.println("Successfully replicated unmodified primer: " + unmodified.seq.equals(expectUnmodified));
        System.out.println("Successfully replicated modified primer: " + modified.seq.equals(expectModified));

    }

}
