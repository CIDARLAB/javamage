package test.Unit;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.DSDNA;
import mage.Tools.FASTA;
import mage.Tools.OutputTools;
import mage.Tools.Pcr.PCR;
import test.Constants;

public class TestOutputTools {

    //used for DSDNA testing
    static String gfp = "atgagcaaaggcgaagaactgtttaccggcgtggtgccgattctggtggaactggatggcgatgtgaacggccataaatttagcgtgagcggcgaaggcgaaggcgatgcgacctatggcaaactgaccctgaaatttatttgcaccaccggcaaactgccggtgccgtggccgaccctggtgaccacctttagctatggcgtgcagtgctttagccgctatccggatcatatgaaacagcatgatttttttaaaagcgcgatgccggaaggctatgtgcaggaacgcaccattttttttaaagatgatggcaactataaaacccgcgcggaagtgaaatttgaaggcgataccctggtgaaccgcattgaactgaaaggcattgattttaaagaagatggcaacattctgggccataaactggaatataactataacagccataacgtgtatattatggcggataaacagaaaaacggcattaaagtgaactttaaaattcgccataacattgaagatggcagcgtgcagctggcggatcattatcagcagaacaccccgattggcgatggcccggtgctgctgccggataaccattatctgagcacccagagcgcgctgagcaaagatccgaacgaaaaacgcgatcatatggtgctgctggaatttgtgaccgcggcgggcattacccatggcatggatgaactgtataaa";

    public static void main(String[] args) throws Exception {
        testGenerateMASCPCRFile();
		//testWritePrimerListsToFile();
        //testGenerateDSDNAPrimerFile();
        //testWriteDSDNAPrimersToFile();
        //testWriteDSDNAPrimersToFileWithLengthError();
        //testGenerateDiversityTableFile();
        //testGetMASCPCRPrimerFileContents();
    }

    public static void testGenerateMASCPCRFile() throws Exception {
        String outfile = Constants.ioDirectory + "/testGenerateMASCPRFile.txt";
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

        System.out.println("generating primers");
        ArrayList<ArrayList<Primer>> primers = pcr.generateAllPrimers(pool);
//		List<List<String>> primers = PCR.getMASCPCRPrimers(pool);

        System.out.println("writing file");
        OutputTools.writePrimersToFile(primers, outfile);
        System.out.println("wrote file " + outfile);
    }

    public static void testWritePrimerListsToFile() throws Exception {
        String outfile = Constants.ioDirectory + "testWritePrimerListsToFile.txt";
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        Oligo insert = Oligo.InsertionFactory(genome, "TTTTTAAAAA", 100, 2, true, "insertion");
        Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730611, 4730612, 2, true, "TAG->TAA");
        Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "deletion");

        List<Oligo> pool = new ArrayList<Oligo>();
        pool.add(insert);
        pool.add(mismatch);
        pool.add(delete);

        System.out.println("generating primers");
		//List<List<String>> primers = PCR.getMASCPCRPrimers(pool);

        System.out.println("writing file");
        List<String> names = new ArrayList<String>();
        names.add("insert");
        names.add("mismatch");
        names.add("delete");

        //OutputTools.writePrimersToFile(names, primers, outfile);
        System.out.println("wrote file " + outfile);
    }

    public static void testGenerateDSDNAPrimerFile() throws IOException {
        String outfile = Constants.ioDirectory + "testGenerateDSDNAPrimerFile.txt";
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        OutputTools.generateDSDNAPrimerFile(genome, gfp, 1000, 2000, outfile);

    }

    public static void testWriteDSDNAPrimersToFileWithLengthError() throws IOException {
        String outfile = Constants.ioDirectory + "testWriteDSDNAPrimersToFileWithLengthError.txt";
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        List<String> primers = DSDNA.getDSDNAPrimers(genome, "", 1000, 2000);
        OutputTools.writeDSDNAPrimersToFile(primers, outfile);
    }

    public static void testWriteDSDNAPrimersToFile() throws IOException {
        String outfile = Constants.ioDirectory + "testWriteDSDNAPrimersToFile.txt";
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        List<String> primers = DSDNA.getDSDNAPrimers(genome, gfp, 1000, 2000);
        OutputTools.writeDSDNAPrimersToFile(primers, outfile);
    }

    public static void testGenerateDiversityTableFile() throws Exception {
        String outfile = Constants.ioDirectory + "testGenerateDiversityTableFile.txt";
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        ArrayList<Oligo> pool = new ArrayList<Oligo>();

        pool.add(Oligo.InsertionFactory(genome, "GCCGCTTTCGCTGTATCCCT", 190, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "AGACAGTCAACAGTAAG", 458, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "GGCCGGA", 2349, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "GTCGATAAGCT", 3599, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "GCTAGAGGAGCGATACGGGATTTAGGAT", 5658, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "GACGATATGAAC", 7900, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "GACTATATA", 14029, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "ATTTGACCT", 15426, 2, true, "oligo"));
        pool.add(Oligo.InsertionFactory(genome, "ATAGCTTTAGGAACCAGACAATGC", 827592, 2, true, "oligo"));

        OutputTools.generateDiversityTrendTableFile(pool, 10, outfile);

    }

    public static void testGetMASCPCRPrimerFileContents() throws Exception {
        String genome = FASTA.readFFN(Constants.blastdirectory, "ecoli.ffn");
        Oligo insert = Oligo.InsertionFactory(genome, "TTTTTAAAAA", 100, 2, true, "insertion");
        Oligo mismatch = Oligo.MismatchFactory(genome, "A", 4730611, 4730612, 2, true, "TAG->TAA");
        Oligo delete = Oligo.DeletionFactory(genome, 89, 100, 2, "deletion");

        List<Oligo> pool = new ArrayList<Oligo>();
        pool.add(insert);
        pool.add(mismatch);
        pool.add(delete);

        //System.out.println(OutputTools.getMASCPCRPrimerFileContents(pool));		
    }
}
