package mage.Switches;

import mage.Tools.BLAST.BlastResult;

public class Blast {

	public static int method = 1;
	
	public static Double score(BlastResult br){
		Double score = 0.0;
		switch (Blast.method) {
			case 1: score = br.bitscore * Math.exp(-1.0 * br.evalue) ; break;
			case 2: score = br.qEnd-br.qStart + 1.0; break;
			default: System.err.println("[Switches] Invalid Blast Scoring System Selected") ; break;
		}
		return score;
	}
	
}
