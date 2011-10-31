package test;

import java.util.ArrayList;

import mage.Oligo;
import mage.Switches;
import tools.Constants;
import tools.FASTA;

public class freeEnergyValidate {

	public static void main (String [] args) throws Exception {
		
		String genome = FASTA.readFFN(Constants.blastdirectory,Oligo.Genome);
		Oligo.buffer_5prime= 35;
		Oligo.buffer_3prime= 25;
		Switches.setFreeEnergyScoringMethod(2);
		System.out.println(genome.length());

		ArrayList <Oligo> pool =  new ArrayList<Oligo>();

		// Test pool of four oligos
//		pool.add(Oligo.InsertionFactory(genome, "ATCG", 190) );
		pool.add(Oligo.InsertionFactory(genome, "TTGG", 458) );
//		pool.add(Oligo.InsertionFactory(genome, "ATCGGCTCGAG", 1408) );
		pool.add(Oligo.InsertionFactory(genome, "GGAATTAACCAA", 2349) );
		
		for (Oligo ol : pool) {
			ol.calc_bg();
			ol.calc_dg();
			System.out.println(ol.getDGasString());
		}
//		ArrayList<String> list = new ArrayList<String>(4);
//		list.add(pool.get(0).getOligo(1));
//		for(String ss: list) {System.out.println(ss);}
//		
//		MFOLD mf = new MFOLD(list);
//		ArrayList<Double> raw_dg = mf.run();
//		
//		for(Double dd : raw_dg) {System.out.println(dd);}
//		
		
		
	}
}