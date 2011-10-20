import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.biojava3.core.sequence.DNASequence;

import tools.BLAST;


/**
 * Oligo object
 * 
 * Similar to a simple Base Pair DNA Sequence with some other tweaks
 * 
 * @author Samir Ahmed
 *
 */
public class Oligo extends DNASequence {

	public static Integer ideal_length = 90;
	public static Integer min_length = 90;

	public static int oligo_count = 0;
	public static int buffer_3prime = 15;
	public static int buffer_5prime = 15;
	
	
	private static HashMap<FeatureIndex,Mistarget> index = new HashMap<FeatureIndex,Mistarget>();

	private String target; 
	private int target_start;
	private int target_end;
	private int target_length;
	
	private int span;
	private int margin;
	private int oligo_min;
	private int oligo_max;
	
	private int feature_count;
	private int oligo_id;
	
	private int genome_start;
	private int genome_end;	

	private ArrayList<Double>  bg_scores;
	private ArrayList<Double>  dg_scores;
	
	public Oligo(String preSequence, String targetSequence, String postSequence, int genomeStart, int genomeEnd){
		super(preSequence+targetSequence+postSequence);
		
		// Store the genome start and end values
		this.genome_start =  genomeStart;
		this.genome_end =  genomeEnd;
		
		// Calculate the span of the oligo
		this.span =  super.getLength();

		// Store the target sequence and length
		this.target = targetSequence;
		this.target_length = targetSequence.length();
		this.target_start = preSequence.length();
		this.target_end =  preSequence.length() + this.target_length;

		// The margin for variation is given by L - 5'Buffer - 3-Buffer - t
		this.margin = Oligo.ideal_length - Oligo.buffer_3prime - Oligo.buffer_5prime - this.target_length;
		
		// Get the index number for the first possible oligo.
		this.oligo_min = 1;
		this.oligo_max = this.oligo_min + this.margin;

		// Set up data members/ structures for keeping track of features.
		this.feature_count = 0;
		this.oligo_id =  ++Oligo.oligo_count;
		
		// Create an ArrayList of Fixed Length for storing blast genome and oligo scores
		this.bg_scores = new ArrayList<Double>(margin); 
		this.dg_scores  = new ArrayList<Double>(margin); 
	}

	
	
	public void calc_bg(){
		try{
			BLAST blast = new BLAST("/Users/mockingbird/dropbox/research/optimization/blast/","genome.ffn");
			HashMap<Integer,String> queries = new HashMap<Integer,String>(); 
			for ( int ii = Oligo.buffer_3prime+1 ; ii < target_start ;ii ++) {
				queries.put(ii,this.getOligo(ii));
				System.out.println(this.getOligo(ii));
			}
			blast.setQuery(queries);
			blast.run();
		}
		catch (Exception ee){ee.printStackTrace();}
	}

	public int getGenomeStart(){
		return this.genome_start;
	}
	
	public int getGenomeEnd(){
		return this.genome_end;
	}
	
	public void addFeature(Mistarget mistarget) {
		super.addFeature(mistarget);
		this.feature_count++;
		index.put(mistarget.getFeatureIndex(), mistarget);
	}
	
	public List<Double> dgList() {
		return this.dg_scores;
	}
	
	public Double dg(int position){
		// If we are in a good range get the following values.
		if (position > Oligo.buffer_5prime && position <= Oligo.buffer_5prime+margin);
		return this.dg_scores.get(position);
	}
	
	public List<Double> bgList() {
		return this.bg_scores;
	}

	public Double bg(int position){
		return this.bg_scores.get(position);
	}
	
	public String getOligo(int start_position, Boolean three_prime_end) throws Exception{

		// If we are inside the 5` or 3` buffer then we throw an error. If we with the target zone, we also through an error
		if ((!three_prime_end) && (start_position <= Oligo.buffer_5prime) || 
			(!three_prime_end) && (start_position > this.target_start) ||
			(three_prime_end) && (start_position <= Oligo.buffer_3prime) || 
			(three_prime_end) && (start_position > this.target_end) ) {
			throw new Exception("Oligo Reference Position Is Too Close to Either the 3 or 5 prime end");
		}
		
		int start;
		
		// If we the position is in reference to the three prime, set the start index according
		if (three_prime_end){
			start = span - start_position - Oligo.ideal_length - Oligo.buffer_3prime; 
		}
		else { start = start_position - Oligo.buffer_5prime; }
		
		// Return the subsequence as a string
		return super.getSubSequence( start  , start +ideal_length).getSequenceAsString();	
	}
	
	public String getOligo(int position) throws Exception {
		return this.getOligo(position, false);
	}
	
	public static Mistarget getMistarget(FeatureIndex fi) throws Exception {
		Mistarget mt;
		if (index.containsKey(fi)){
			mt = index.get(fi);
		}
		else
		{
			throw new Exception("No such mismatch");
		}
		return mt;
	}
	
	

	/**
	 * To get a new index for a feature - the getNewFeatureIndex method will return a FeatureIndex
	 * That may be used to create a new feature.
	 * 
	 * @return
	 */
	public FeatureIndex getNewFeatureIndex(){
		this.feature_count++;
		return new FeatureIndex(this.oligo_id,this.feature_count);
	}

	//	public Mistarget getMistargetByIndex(int index_no){
	//		for ( FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> mt: this.getFeatures() ){
	//			
	//		}
	//	}
	//	
	public int getFeatureCount(){
		return this.feature_count;
	}

	public int getOligoId(){
		return this.oligo_id;
	}

	public static void main(String[] args) {
		Oligo ol= new Oligo("ATTGCATCGATACATAAGATGTCTCGACCGCATGCGCAACTTGTGAAGTGTCTACTATCC",
				"GGGATTATTATTGGG","CTAAGCCCATTTCTCGCACAATAACCCCTGAATGTGTCCGCATCTGATGTTACCCGGGTT",1,100);
		
		ol.calc_bg();
	}

}
