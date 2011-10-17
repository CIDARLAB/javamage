import org.biojava3.core.sequence.DNASequence;


public class test {
	
	public static void main (String[] args){
		
		DNASequence seq = new DNASequence("atcgATCGCGATGACGACGAGGACGACAGCAGACGACGACGAGTC");
		seq.addFeature( new Mistarget(0,5,4.5));
		
		System.out.println(seq.getFeatures().get(0).getLocations().getEnd());
		
	}
}
