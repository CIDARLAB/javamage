package mage;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

import mage.Switches.OligoScore;

import org.biojava3.core.sequence.DNASequence;

import tools.BLAST;
import tools.BLAST.BlastResult;
import tools.Constants;
import tools.FASTA;
import tools.MFOLD;


/**
 * Oligo object
 * 
 * Similar to a simple Base Pair DNA Sequence with some other tweaks
 * 
 * @author Samir Ahmed
 *
 */
public class Oligo extends DNASequence {


	public	static String Genome = "genome2.ffn";
	public 	static Integer ideal_length = 90;
	public 	static Integer min_length = 90;

	public 	static int oligo_count = 0;
	public 	static int buffer_3prime = 15;
	public 	static int buffer_5prime = 15;
	public 	static ArrayList<Oligo>	all = new ArrayList<Oligo>();

	public 	static HashMap<Integer,Oligo> oligo_map = new HashMap<Integer,Oligo>();
	private static HashMap<FeatureIndex,Mistarget> index = new HashMap<FeatureIndex,Mistarget>();

	final private String target; 
	final private int target_start;
	final private int target_end;
	final private int target_length;

	final private int span;
	final private int margin;
	final private int oligo_min;
	final private int oligo_max;

	private int feature_count;
	final private int oligo_id;

	private int primary_position;

	private int x;
	final private int genome_start;
	final private int genome_end;	

	private ArrayList<Double>  	bg_scores;
	private ArrayList<Double>  	dg_scores;
	private ArrayList<Double>  	dg_raw_scores;
	private ArrayList<Double>  	bo_scores;
	private ArrayList<Integer>	valid_dg_positions;
	private ArrayList<Integer>	bo_sorted;
	private ArrayList<Integer>	bg_sorted;
	private ArrayList<Integer>	dg_sorted;
	private ArrayList< LinkedList<Mistarget> > bg_mistargets;

	// Set of Associated Mistargets with the given span
	public ArrayList< Mistarget > mt_collection;

	// Set of Valid Mistargets with the given subsequence of the span
	public ArrayList< Mistarget> valid_mt;

	// Datamembers that indicate the optimized oligo
	private String 	optimized;
	private int		opt_start;
	private	int		opt_end;
	private Double bo_weighted;
	
	// Scoring Objects
//	private OligoScore	optimizedScore;
	private OligoScore 	currentScore;
	

	/**
	 *	This is a static Factory Method that creates an Oligo for Insertion Mutations
	 * 
	 * @param genome			A string containing the genome to be modified by MAGE
	 * @param target			A string containing the target mutation to be inserted/mismatched
	 * @param targetPosition	Integer containing the starting position
	 * 
	 * @return
	 * @throws Exception 
	 */
	public static Oligo InsertionFactory(String genome, String target, int targetPosition) throws Exception {

		if (genome.length() > (targetPosition+Oligo.ideal_length-Oligo.buffer_3prime-1 - target.length())  ) {

			// Define the Genome Starting Position
			int genome_start = targetPosition+Oligo.buffer_3prime - Oligo.ideal_length + target.length() + 1 ;

			// Extract a Subsequence from that point - THIS FOR A STRING i.e INDEXED FROM ZERO
			String preSequence = genome.substring( genome_start -1 ,targetPosition);

			// Define a Genome ending position
			int genome_end = targetPosition +(Oligo.ideal_length-Oligo.buffer_5prime-target.length());

			// Extract a Subsequence from the target Position to the end. THIS IS ALSO FOR A STRING i.e INDEXED FROM ZERO
			String postSequence = genome.substring(targetPosition , genome_end);

			// Return the new Oligo that was just made
			return new Oligo(preSequence, target, postSequence, genome_start, genome_end);	

		}
		else {
			//return null;
			System.out.println("Genome Size = " +genome.length() +  "Requires " +targetPosition + " + " + (Oligo.ideal_length - Oligo.buffer_3prime) ) ;
			throw new Exception("[Insertion Factory] Oligo not defined on correct range" );			
		}

	}

	public Oligo(String preSequence,  String targetSequence, String postSequence, int genome_start, int genome_end) throws Exception{
		super(preSequence+targetSequence+postSequence);

		// Store the genome start and end values
		this.genome_start 	= genome_start;
		this.genome_end 	= genome_end;

		// Calculate the span of the oligo
		this.span 			= super.getLength();
		this.x 				= preSequence.length()+1;

		// Store the target sequence and length
		this.target 		= targetSequence;
		this.target_length 	= targetSequence.length();
		this.target_start	= preSequence.length()+1;
		this.target_end 	= preSequence.length() + this.target_length;

		// The margin for variation is given by L - 5'Buffer - 3-Buffer - t
		this.margin 		= Oligo.ideal_length - Oligo.buffer_3prime - Oligo.buffer_5prime - this.target_length +1;

		// Get the index number for the first possible oligo.
		this.oligo_min 		= 1;
		this.oligo_max 		= this.oligo_min + this.margin;

		// Set up data members/ structures for keeping track of features.
		this.feature_count 	= 0;
		this.oligo_id 		= ++Oligo.oligo_count;

		// Create arrayList of that hold a list of all the mis-targets in a blast genome sequence
		this.bg_mistargets 	= new ArrayList< LinkedList<Mistarget>>(this.margin); 

		// Create an ArrayList of Fixed Length for storing blast genome and oligo scores
		this.bg_scores 			= new ArrayList<Double>(this.margin);
		this.dg_scores  		= new ArrayList<Double>(this.margin);
		this.dg_raw_scores  	= new ArrayList<Double>(this.margin);
		this.bo_scores			= new ArrayList<Double>(this.margin);
		this.valid_dg_positions = new ArrayList<Integer>(this.margin);
		this.bo_sorted			= new ArrayList<Integer>(this.margin);
		this.bg_sorted			= new ArrayList<Integer>(this.margin);
		this.dg_sorted			= new ArrayList<Integer>(this.margin);
		
		// Set the primary position to -1 until it has been determined
		this.primary_position = -1;

		// Create empty Mistarget Sets
		this.mt_collection			= new ArrayList <Mistarget> ();
		this.valid_mt		= new ArrayList <Mistarget> ();

		// Set the score of the weighted BO to zero
		this.bo_weighted	= 0.0;

		// Add this span to the oligo Map
		Oligo.oligo_map.put(this.oligo_id, this);

		// Take the first oligo as the optimized Oligo
		this.optimized 		= getOligo(this.oligo_min);
		this.calcOptimizedBounds(this.oligo_min);

		// Add the oligo to the collection of all oligos
		Oligo.all.add(this);
		
		// Create scoring values;
//		this.optimizedScore 	= new OligoScore();
		this.currentScore		= new OligoScore();
		
		// Notify when Oligo Is Recorded
		System.out.println("Oligo ID: " + this.oligo_id + "; Span : "+this.span+"; Margin: "+this.margin+"; Target : " + this.target );
	}

	public void setOptimized(int start_position) throws Exception{
		this.optimized = getOligo(start_position);
		calcOptimizedBounds(start_position);
	}
	
	private void calcOptimizedBounds(int start_position){
		this.opt_start	= start_position;
		this.opt_end	= start_position+Oligo.ideal_length-1;
	}

	/**
	 * 	Calculates the primary score associated with the span.  Use the getPrimaryScoreAsString method to retreive it
	 */
	public void calc_primaryScore() {
		this.primary_position = mage.Switches.PrimaryScore.score(bg_scores,dg_scores) + this.oligo_min;	
	}

	/**
	 * The primary score is the optimal free energy and blast genome score combined
	 * @return	A formatted string in teh form "Oligo_%d Primary Score ( %.3f , %.3f )"
	 */
	public String getPrimaryScoreAsString(){
		return String.format("Oligo_%d Primary Score ( %.3f , %.3f )",this.oligo_id , this.bg(this.primary_position),this.dg(this.primary_position)); 
	}

	/**
	 * Calculates the Free energy for all positions on the span
	 * 
	 * Stores result in the dg_scores arraylist<Double>
	 * 
	 */
	public void calc_dg(){
		try{
			ArrayList <String> list = new ArrayList<String>(this.margin);
			for ( int ii = this.oligo_min ; ii < this.oligo_max; ii++){
				list.add(this.getOligo(ii));
			}
			
			// Create an new instance of MFOLD using the list of oligos
			// Run Mfold and capture results in the dg_raw_scores list
			MFOLD mfold = new MFOLD(list);
			dg_raw_scores = mfold.run();			

			// Calculate the score with the switch then add it to the valid dg positions list
			// NOTE: That valid dg positions starts indexing from 1 not zero
			int counter = 1;
			for (Double dg_value : dg_raw_scores){ 
				Double score = mage.Switches.FreeEnergy.score(dg_value);
				dg_scores.add(score);
				
				// Check if we fall below the threshold, if so we add this to the valid dg positions list
				if (mage.Switches.FreeEnergy.threshold(score)) {
					valid_dg_positions.add(counter);
				}
				dg_sorted.add(counter);
				counter++;
			}
			
			Collections.sort(dg_sorted, new Comparator<Integer>() {
				//Implementing stanard compare fucntion for sorting ascending
				public int compare(final Integer ii, final Integer jj) {
					return (int) ( (dg_scores.get(ii-1))- (dg_scores.get(jj-1)));
				}
			});

			
		}
		catch (Exception ee) {ee.printStackTrace();}
	}


	/**
	 * 
	 * Calculates the Blast genome Score for all the oligos in the span
	 * 
	 * Stores the result in the bg_score arrayList<Double>
	 * 
	 */
	public void calc_bg(){
		try{
			BLAST blast = new BLAST(Constants.blastdirectory,Oligo.Genome);
			HashMap<Integer,String> queries = new HashMap<Integer,String>(); 

			ArrayList< ArrayList<Double>> score_list = new ArrayList< ArrayList<Double> > (this.margin) ;

			// Starting from leftmost position on the span, move to the right and get
			for ( int ii = this.oligo_min ; ii < this.oligo_max ; ii ++) {

				// Extract an oligo for every posible position and put it in the map
				queries.put(ii,this.getOligo(ii));
				System.err.println(this.getOligo(ii));
				score_list.add(new ArrayList<Double>());
			}

			/* Set the BLSAT Query and then RUN IT, capturing the results inside a LIST*/
			blast.setQuery(queries);
			List<BlastResult> br = blast.run();

			for ( BLAST.BlastResult result :  br) {

				/* Check if we have a mistarget not targeted alignment*/
				if( !( (result.sStart >= this.genome_start) && (result.sStart <= this.genome_end) ) ||
						!( (result.sEnd >= this.genome_start) && (result.sEnd <= this.genome_end) ) 	){

					/* Calculate the score using for the BLAST Genome : SWITCH */
					Double score = mage.Switches.Blast.score(result);

					/* Store the score */
					score_list.get(result.oligoID - (this.oligo_min)).add(score);

					//					FeatureIndex = new FeatureIndex();
					//					Mistarget mt = new Mistarget( , , ,score));
				}
			}

			int count = 0;

			// If we have an empty array we simply put in the score zero, otherwise we put in the Genome Score per position
			for (ArrayList<Double> position : score_list){
				if (position.isEmpty()){ bg_scores.add( count ,0.0); }
				else { bg_scores.add( count, mage.Switches.BlastGenomeTotal.score(position)); }
				count ++;
				bg_sorted.add(count);
			}
			
			// Create a list of sorted Indices based on bg Scores
			Collections.sort(bg_sorted, new Comparator<Integer>() {
				//Implementing stanard compare fucntion for sorting ascending
				public int compare(final Integer ii, final Integer jj) {
					return (int) ( (bg_scores.get(ii-1))- (bg_scores.get(jj-1)));
				}
			});

			
			
		}
		catch (Exception ee){ee.printStackTrace();}
	}

	/**
	 * Returns the start position of the oligo that represents the primary position
	 * @return	Integer containing the starting position of the primary oligo.
	 */
	public int getPrimaryPosition() {
		return this.primary_position;
	}

	public String getBGasString(){	
		StringBuilder sb =  new StringBuilder();

		for (int ii = this.oligo_min, jj=0 ; ii< oligo_max ; ii++, jj++) {
			sb.append(String.format( "%d \t \t %.3f \n",jj, bg(ii)) ); 
		}

		return sb.toString();
	}

	public String getDGasString(){

		StringBuilder sb =  new StringBuilder();
		for (int ii = this.oligo_min, jj=0 ; ii< oligo_max ; ii++, jj++) {
			sb.append(String.format( "%d \t \t %.3f \n",jj, dg(ii) )); //, dg_raw_scores.get(ii-this.oligo_min) )); 
		}

		return sb.toString();
	}

	public int getGenomeStart(){
		return this.genome_start;
	}

	public List<Double> dgList() {
		return this.dg_scores;
	}

	private Double dg(int position){
		if (position >= this.oligo_max && position < this.oligo_min);
		return this.dg_scores.get(position-this.oligo_min);
	}

	public List<Double> bgList() {
		return this.bg_scores;
	}

	public String getTarget(){
		return this.target;
	}

	public int getGenomeEnd(){
		return this.genome_end;
	}

	private Double bg(int position){
		return this.bg_scores.get(position-this.oligo_min);
	}

	/**
	 * Returns an oligo that is the subset of the span.
	 * 
	 * 
	 * e.g <pre> 
	 * {@code
	 * oligo.getOligo(1)	//will return an oligo of the first n base pairs for n = ideal length of oligo. 
	 * oligo.getOligo(0)	//will throw an error, because the first position is 1.
	 * oligo.getOligo(314)	//will throw an error, given that it indexes outside the range of the span of the oligo. 
	 * }
	 * </pre>
	 * @param start_position	An integer representing the starting position of the oligo of ideal length
	 * @return					An oligo of ideal length
	 * @throws Exception		If the oligo is outside the bounds of the span, and error is thrown
	 */
	private String getOligo(int start_position) throws Exception{

		// Check if this an oligo with the given starting position can be extracted
		if (start_position > this.oligo_max || start_position < this.oligo_min){
			throw new Exception("Extracting Oligo outside of Span");
		}

		// Return the subsequence as a string
		return super.getSubSequence( start_position, start_position +ideal_length-1).getSequenceAsString();	
	}


	public void addMistarget(Mistarget mt) {
		this.mt_collection.add(mt);
	}


	/** 
	 * Shift Optimized takes the current span, extracts a subsequence from the desired starting position
	 * <p>
	 * It creates and oligo of ideal length and then updates the mistargets associated with the oligo
	 * </p>
	 * 
	 * @param start_position
	 * @throws Exception
	 */
	public void setOligo(int new_start_position) throws Exception {
		// Get new optimized Oligo
		this.setOptimized(new_start_position);

		this.valid_mt.removeAll(this.valid_mt);
		// For each associated mistarget check if it is valid
		for (Mistarget mt : mt_collection){
			if (mt.isValid(this)) {
				this.valid_mt.add(mt);
			}
		}
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

	/**
	 * Returns the number of features associated with this oligo
	 * @return	Returns an positive integer representing the number of linked mistarget
	 */
	public int getMistarget(){
		return this.valid_mt.size();
	}

	/**
	 * Every oligo has a unique static id number
	 * @return	An integer holding the unique id number of the referenced oligo
	 */
	public int getOligoId(){
		return this.oligo_id;
	}

	/**
	 * Get the optimized Oligos end position on the span
	 * @return	End position relative the start of the oligo's span
	 */
	public int getOptimizedEnd() {
		return this.opt_end;
	}

	/**
	 * Gets the optimized Oligos start position on the span
	 * @return	Start Position relative the start of the oligo's span
	 */
	public int getOptimizedStart() {
		return this.opt_start;
	}


	public double getWeightedBOScore() {
		return this.bo_weighted;
	}

	/**
	 * Blast two oligos against each other. Both oligos will have the associated mistargets registerd.
	 * 
	 * @param olA
	 * @param olB
	 * @throws IOException
	 */
	public static void BlastOligo(Oligo subject, List<Oligo> queries) throws IOException {

		// Create a .FFN with from the Subject
		String filename = "oligo.ffn";
		String directory = Constants.bodirectory;
		FASTA.writeFFN(directory, filename, subject.getSequenceAsString().toString());

		// Create a BLASTDB with .FFN files
		BLAST bl = new BLAST(directory,filename);

		// Make a map queries
		HashMap<Integer,String> query= new HashMap<Integer,String>();
		for (Oligo ol : queries) { query.put(ol.getOligoId(), ol.toString().toString()) ;}

		// Assign the queries to this instance of blast
		bl.setQuery(query);

		// Run BLAST
		List<BlastResult> results = bl.run();

		// Create a mistarget for every blast result
		for (BlastResult br:results) {
			Mistarget.mistarget_collection.add( new Mistarget(br,subject.getOligoId(),br.oligoID ));
		}

	}

	/**
	 * Calculates the BO Score, weighted for the optimized oligo
	 * 
	 * This uses switches.BlastOligoWeight.score 
	 * @throws Exception 
	 * 
	 */
	public void calc_primary_bo() {

		// Calculate and assign the blast oligo score
		this.bo_weighted= mage.Switches.BlastOligo.score(this);

	}

	/**
	 * Calculates the Blast Oligo scores for every oligo on the genome given the current valid set
	 * @throws Exception
	 */
	public void calc_bo() throws Exception {
		
		// Set every possible oligo on the span as the optimized oligo, then calculate the associated BO Value
		for ( int ii = this.oligo_min; ii < this.oligo_max; ii++ ) {
			this.setOligo(ii);
			Double score = mage.Switches.BlastOligo.score(this);
			this.bo_scores.add(score);
			
			// Add the Index Number for sorting
			this.bo_sorted.add(ii);
			System.err.println(this.optimized);
		}
		
		// Sort Ascending by BO value
		Collections.sort(this.bo_sorted, new Comparator <Integer>() {
			
			//Implementing stanard compare fucntion for sorting ascending
			public int compare(final Integer ii, final Integer jj) {
				return (int) ( (bo_scores.get(ii-1))- (bo_scores.get(jj-1)));
			}
		});
		
		// Reset the Oligo to the primary position
		this.setOligo(this.primary_position);
	}
	
	/**
	 * Sorts a given pool of oligos by WeightedBO Score
	 * 
	 * @param pool	An Array List of Oligos. After calling method, List will be in sorted order
	 */
	public static void sort(ArrayList<Oligo> pool){

		// Now sort the mt_collection by 
		Collections.sort(pool, new Comparator <Oligo>(){

			// Implementing standard compare function by 
			public int compare(final Oligo o1, final Oligo o2) { 
				return (int) (o2.getWeightedBOScore() - o1.getWeightedBOScore()); 
			}
		});
	}

	/**
	 * Returns the score at the current position;
	 * @return	Score in the form of type OligoScore
	 */
	public OligoScore currentScore() {
		
		// Set the current score to the score at the optimized position
		this.currentScore = this.scoreAt(this.opt_start);
		
		// Return the current score
		return this.currentScore;
	}
		
	/**
	 * Returns the OligoScore at the indicated position	
	 * @param start_position
	 * @return
	 */
	public OligoScore scoreAt(int start_position) {
		
		// Create a new Oligo Score and Add the Score Components Individually
		OligoScore newScore = new OligoScore();
		try {
		
			newScore.BlastOligo(	bo_scores.get(start_position-1));
			newScore.BlastGenome(	bg_scores.get(start_position-1));
			newScore.FreeEnergy(	dg_scores.get(start_position-1));
		}
		catch (Exception ee) {
			ee.printStackTrace();
		}
		
		// Return the score
		return newScore;
	}

	/**
	 * Get the first possible starting position on the margin
	 * @return	Integer with starting index
	 */
	public int getOligoMinimum() {return this.oligo_min;}

	/**
	 * Get the last possible starting position on the margin
	 * @return	Integer with starting index
	 */
	public int getOligoMaximum() {return this.oligo_max;}

	/**
	 * Get the number for possible oligos on the span i.e the length of the margin
	 * @return Integer count of possible oligos
	 */
	public int getMarginLength() {return this.margin;}

	
//	/**
//	 * Returns a stack containing a prioritized Integer values;
//	 * @return	Sorted list of integers containing oligo positions
//	 */
//	public Stack<Integer> getBOStack() { 
//		Stack<Integer> stack = new Stack<Integer>();
//		for (Integer ii: this.bo_sorted) {
//			stack.push(ii);
//		}
//		return stack;
//	}

	/**
	 * Get A List of Indicies corresponding to sorted BG Positions
	 * @return
	 */
	public List<Integer> getBGSorted() { return this.bg_sorted; }
	
	/**
	 * Get A Lst of Indicies corresponding to sorted BO positions
	 * @return
	 */
	public List<Integer> getBOSorted() {return this.bo_sorted; }
	
	/**
	 * Gets a List of Indicices corresponding to sorted DG positions
	 * @return
	 */
	public List<Integer> getDGSorted() { return this.dg_sorted; }
	
	/**
	 * Get a list of Indicies corresponding to DG positions below the designated threshold
	 * @return
	 */
	public List<Integer> getValidDG() { return this.valid_dg_positions;}
	
}
