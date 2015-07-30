package mage.Core;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import mage.Tools.BLAST;
import mage.Tools.BLAST.BlastResult;
import mage.Tools.Constants;
import mage.Tools.FASTA;
import mage.Tools.MFOLD;
import mage.Tools.SequenceTools;
import org.biojava3.core.sequence.DNASequence;




/**
 * Oligo object
 * 
 * Similar to a simple Base Pair DNA Sequence with some other tweaks
 * 
 * @author Samir Ahmed
 *
 */
public class Oligo extends DNASequence {


	public 	static ArrayList<Oligo>	all = new ArrayList<Oligo>();
	public 	static int buffer_3prime = 15;	
	public 	static int buffer_5prime = 15;
	public	static String	Directory = Constants.workingdirectory;

	public	static String 	Genome = "genome.ffn";
	public 	static Integer ideal_length = 90;
	public 	static Integer min_length = 90;
	private static int oligo_count = 0;

	public 	static HashMap<Integer,Oligo> oligo_map = new HashMap<Integer,Oligo>();

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
		String directory = Constants.workingdirectory;
		FASTA.writeFFN(directory, filename, subject.getSequenceAsString().toString());

		// Create a BLASTDB with .FFN files
		BLAST bl = new BLAST(directory,filename);

		// Make a map queries
		HashMap<Integer,String> query= new HashMap<Integer,String>();
		for (Oligo ol : queries) { 
			query.put(ol.getOligoId(), SequenceTools.ReverseCompliment(ol.getSequenceAsString())) ;
		}

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
	 * Oligo Factory,
	 * 
	 * Given a genome, and a target.  The oligo factory will return the desired Oligo object
	 * defined by the parameters of the target object.
	 * 
	 * @param genome	A genome of sufficient length
	 * @param tt		A target mutation object 
	 * @return 			An unoptimized oligomer to perform the target mutation
	 * @throws Exception
	 */
	public static Oligo OligoFactory(String genome, Target tt) throws Exception {
		//remove whitespace
		genome = genome.replaceAll("\\s", "");
		String sequence = tt.sequence.replaceAll("\\s","");
		
		if ( tt.type.startsWith("I") ){
			return Oligo.InsertionFactory(genome, sequence, tt.left_position, tt.replichore , tt.sense , tt.gene_name);
		}
		else if (tt.type.startsWith("M")) {
			return Oligo.MismatchFactory(genome, sequence, tt.left_position, tt.right_position, tt.replichore, tt.sense, tt.gene_name) ;
		}
		else if (tt.type.startsWith("D")) {
			return Oligo.DeletionFactory(genome, tt.left_position, tt.right_position, tt.replichore, tt.gene_name) ;
		}
		else {
			System.out.println("Could Not Create Oligo");
			throw new Exception("[OligoFactory] Oligo not defined on correct range" );	
		}

	}

	/**
	 * Defines a mismatch on a surface.
	 * @param genome 		The genome from which the oligo is being made
	 * @param target		The target sequence that will be put onto the genome
	 * @param left_position	1-indexed, inclusive. The left position, indicating the point at which the mutation will be inserted into the genome
	 * @param right_position 1-indexed, exclusive. The end position of where the oligo is to cut off
	 * @param replichore	The replichore 1 or 2
	 * @param sense			True if sense is postive and false if negative
	 * @param name			The name of this target
	 * @return				Returns a oligo defined by the input parameters
	 * 
	 * @throws Exception	In the event that oligo is incompletely/improperly specified
	 */
	
	public static Oligo MismatchFactory(String genome, String target, int left_position, int right_position, int replichore, boolean sense, String name)  throws Exception { 
		if (genome.length() >  (right_position+Oligo.ideal_length - Oligo.buffer_5prime -1) && (left_position > (Oligo.ideal_length - Oligo.buffer_3prime))) 
		{
			// Define the starting position on the genome)
		//correct for 0-indexing
                    left_position--;
                    right_position--;
                    
			int genome_start = left_position - Oligo.ideal_length + target.length() + Oligo.buffer_3prime -1;

			// Pulls from the genome string, this start position -1 up to left position -1 (not inclusive)
			String preSequence =  genome.substring(genome_start, left_position);
                        
			// Define the ending position on the genome
			int genome_end  = right_position + Oligo.ideal_length - target.length() - Oligo.buffer_5prime;

			// Pulls the genome string, the start position is   right_position (inclusive), genome_end (exclusive)
			String postSequence = genome.substring(right_position, genome_end ); 

                        // Take the reverse compliments of genome for -1,+1 ... ignore for -2,+2
			if (replichore == 1){

				// Get RC
				String reverseComp  = mage.Tools.SequenceTools.ReverseCompliment(preSequence+postSequence);

				// Calculate index to split from and then reassign post and pre sequnence
				int splitIndex 		= postSequence.length();
				postSequence 		= reverseComp.substring(splitIndex);
				preSequence 		= reverseComp.substring(0,splitIndex);
			}
			if ( ((replichore==2) && !sense) || ((replichore==1) && sense) ) {
				target = mage.Tools.SequenceTools.ReverseCompliment(target);	
			}
			
			// Return the new Oligo that was just made
			return new Oligo(preSequence, target, postSequence, genome_start, genome_end, OligoType.MISMATCH, name);	

		}
		else {
			//return null;
			System.out.println("Could Not Create Oligo");
			throw new Exception("[Mismatch Factory] Oligo not defined on correct range" );			
		}

	}

	/**
	 * Deletion factory generates and an oligo that represents a deletion
	 * The bases from left_position to right_postion (inclusive) will be removed
	 * For example, DeletionFactory("ABCDEFG", 3, 5,...) will generate "ABEFG"
	 * @param genome 		The genome from which the oligo is being made
	 * @param left_position	1-indexed left position, inclusive. Indicating the point at which the mutation will start
	 * @param replichore	The replichore 1 or 2
	 * @param right_position 1-indexed, exclusive. An Integer to specify the position after which the genome is to resume
	 * @param name			The name of this oligo
	 * @return				Returns a oligo defined by the input parameters
	 * @throws Exception	In the event that oligo is incompletely/improperly specified
	 */
	//TODO: Confirm that the indexing on left and right positions is intuitive.
	//Currently, both coordinates are INCLUSIVE- the change will affect the bases at both indicated indexes
	public static Oligo DeletionFactory(String genome, int left_position, int right_position, int replichore, String name) throws Exception {
	//correct for 0-indexing
                    left_position--;
                    right_position--;
                              
            if (genome.length() >  (right_position+Oligo.ideal_length - Oligo.buffer_5prime) && (left_position > (Oligo.ideal_length - Oligo.buffer_3prime))) {

			// Define the starting position on the genome)
			int genome_start = left_position - Oligo.ideal_length + Oligo.buffer_3prime;

			// Pulls from the genome string, this start position -1 up to left position -1 +1 (not inclusive)
			String preSequence =  genome.substring(genome_start, left_position);

			// Define the ending position on the genome
			int genome_end  = right_position + Oligo.ideal_length - Oligo.buffer_5prime;

			// Pulls the genome string, from right_position (-1 due to 0 indexing, genome_end (non-inclusive)
			String postSequence = genome.substring(right_position, genome_end); 

			// Take the reverse compliments of genome for -1,+1 ... ignore for -2,+2
			if (replichore == 1){

				// Get RC
				String reverseComp  = mage.Tools.SequenceTools.ReverseCompliment(preSequence+postSequence);

				// Calculate index to split from and then reassign post and pre sequnence
				int splitIndex 		= preSequence.length();
				postSequence 		= reverseComp.substring(splitIndex);
				preSequence 		= reverseComp.substring(0,splitIndex);
			}
			return new Oligo(preSequence,"",postSequence,genome_start,genome_end, OligoType.DELETION, name);
		}
		else {
			//return null;
			System.out.println("Could Not Create Oligo") ;
			throw new Exception("[Deletion Factory] Oligo not defined on correct range");			
		}

	}

	/**
	 * Generates an insertion oligo given a set of parameters
	 * 
	 * @param genome 		The genome from which the oligo is being made
	 * @param target		The target sequence that will be put onto the genome
	 * @param left_position	The left position (1-indexed,inclusive), indicating the point at which the mutation will be inserted into the genome
	 * @param replichore	The replichore 1 or 2
	 * @param sense			True if sense is postive and false if negative
	 * @param name			The name of this target
	 * @return				Returns a oligo defined by the input parameters
	 * 
	 * @throws Exception	In the event that oligo is incompletely/improperly specified
	 */
	public static Oligo InsertionFactory(String genome, String target, int left_position, int replichore, boolean sense, String name) throws Exception {
		//correct for 0-indexing
                    left_position--;	
            if ( (genome.length() > (left_position+Oligo.ideal_length-Oligo.buffer_3prime - target.length())  ) && (left_position > 60 ) ) {

			// Define the starting pisition on the genome)
			int genome_start = left_position + target.length() - Oligo.ideal_length + Oligo.buffer_3prime;

			// Pulls from the genome string, this start position -1 up to left position -1 (not inclusive)
			String preSequence =  genome.substring(genome_start, left_position);

			// Define the ending position on the gneome
			int genome_end  = left_position + Oligo.ideal_length - Oligo.buffer_5prime -target.length();

			// Extract a Subsequence from the target Position to the end. THIS IS ALSO FOR A STRING i.e INDEXED FROM ZERO
			String postSequence = genome.substring(left_position , genome_end);

			// Take the reverse compliments of genome for -1,+1 ... ignore for -2,+2
			if (replichore == 1){

				// Get RC
				String reverseComp  = mage.Tools.SequenceTools.ReverseCompliment(preSequence+postSequence);

				// Calculate index to split from and then reassign post and pre sequnence
				int splitIndex 		= postSequence.length();
				postSequence 		= reverseComp.substring(splitIndex);
				preSequence 		= reverseComp.substring(0,splitIndex);
			}

			if ( ((replichore==2) && !sense) || ((replichore==1) && sense) ) {
				target = mage.Tools.SequenceTools.ReverseCompliment(target);	
			}
			// Return the new Oligo that was just made
			return new Oligo(preSequence, target, postSequence, genome_start, genome_end, OligoType.INSERTION,name);	


		}
		else {
			//return null;
			//TODO: handle the case where insertion is near the start of the sequence
			System.out.println("Could Not Create Oligo");
			throw new Exception("[Insertion Factory] Oligo not defined on correct range" );			
		}

	}

	/**
	 * Sorts a given pool of oligos by WeightedBO Score
	 * 
	 * @param pool  An Array List of Oligos. After calling method, List will be in sorted order
	 */
	public static void sort(List<mage.Core.Oligo> pool){


		// Now sort th by 
		Collections.sort(pool, new Comparator<Oligo>(){

			// Implementing standard compare function by 
			public int compare(final Oligo o1, final mage.Core.Oligo o2) { 

				double o1_min= o1.getGreedyScore();
				double o2_min= o2.getGreedyScore();

				return (int) (o2_min - o1_min);
			}
		});
	}


	private ArrayList<Double>  	bg_scores;
	private ArrayList<Integer>	bg_sorted;
	private ArrayList<Double>  	bo_scores;

	private ArrayList<Integer>	bo_sorted;
	private Double 				bo_weighted;
	private OligoScore 			currentScore;

	private ArrayList<Double>  	dg_raw_scores;
	private ArrayList<Double>  	dg_scores;
	private ArrayList<Integer>	dg_sorted;

	final 	private int genome_end;
	final 	private int genome_start;
	private int			greedy_choice;
	final 	private int margin;

	// Set of Associated Mistargets with the given span
	public 	ArrayList< Mistarget > 	mt_collection;
	private ArrayList<Integer>		valid_dg_positions;
	public 	ArrayList< Mistarget> 	valid_mt;

	final 	private	int		oligo_id;
	public final 	 int		oligo_max;
	public final 	 int		oligo_min;

	private	int				opt_end;
	private int				opt_start;
	private String 			optimized;
	private int 			primary_position;
	private int				optMagePosition;
	
	//is the oligo used for insertion, deletion, or mismatch?
	private OligoType oligoType;

	// Immutable members
	final 	public String 	sequence;
	final 	public int 	span;
	final 	public String 	target;
	final 	public int 	target_length;
	final   public int target_position;

	final public String 	name;

	//public final int x; //not used, and incorrectly calculated if the oligo is reversed
	private ArrayList<OligoScore> scores;

	/**
	 * Makes an oligo defined by the genome starting and ending points
	 * 
	 * @param preSequence
	 * @param targetSequence
	 * @param postSequence
	 * @param genome_start
	 * @param genome_end
	 * @param name
	 * @throws Exception
	 */
	public Oligo(String preSequence,  String targetSequence, String postSequence, int genome_start, int genome_end, OligoType oligoType, String name) throws Exception{
		super(preSequence+targetSequence+postSequence);

		this.name = name.replaceAll("\\s+", "");

		this.sequence = preSequence+targetSequence+postSequence;
		// Store the genome start and end values
		this.genome_start 	= genome_start;
		this.genome_end 	= genome_end;

		// Calculate the span of the oligo
		this.span 			= super.getLength();

		//this.x 				= preSequence.length()+genome_start;

		// Store the target sequence and length
		this.target 		= targetSequence;
		this.target_length 	= targetSequence.length();
		this.target_position= preSequence.length();
		// The margin for variation is given by L - 5'Buffer - 3-Buffer - t
		this.margin 		= Oligo.ideal_length - Oligo.buffer_3prime - Oligo.buffer_5prime - this.target_length +1;
		
		// Get the index number for the first possible oligo.
		this.oligo_min 		= 1;
		this.oligo_max 		= this.oligo_min + this.margin;

		// Set up data members/ structures for keeping track of features.
		this.oligo_id 		= ++Oligo.oligo_count;

		// Create an ArrayList of Fixed Length for storing blast genome and oligo scores
		this.bg_scores 			= new ArrayList<Double>(this.margin);
		this.dg_scores  		= new ArrayList<Double>(this.margin);
		this.dg_raw_scores  	= new ArrayList<Double>(this.margin);
		this.bo_scores			= new ArrayList<Double>(this.margin);
		this.valid_dg_positions = new ArrayList<Integer>(this.margin);
		this.bo_sorted			= new ArrayList<Integer>(this.margin);
		this.bg_sorted			= new ArrayList<Integer>(this.margin);
		this.dg_sorted			= new ArrayList<Integer>(this.margin);
		this.scores				= new ArrayList<OligoScore>(this.margin);
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
		this.greedy_choice= this.oligo_min;
		this.optMagePosition = -1;
		//this.calcOptimizedBounds(this.oligo_min);

		//set the OligoType to allow replacement efficiency calculation
		this.oligoType = oligoType;
		
		// Add the oligo to the collection of all oligos
		Oligo.all.add(this);

		// Create scoring values;
		//		this.optimizedScore 	= new OligoScore();
		this.currentScore		= new OligoScore();

		// Notify when Oligo Is Recorded
		//System.out.println("Oligo ID: " + this.oligo_id + "; Span : "+this.span+"; Margin: "+this.margin+"; Target : " + this.target );
	}

	public void addMistarget(Mistarget mt) {
		this.mt_collection.add(mt);
	}

	private Double bg(int position){
		return this.bg_scores.get(position-this.oligo_min);
	}

	public List<Double> bgList() {
		return this.bg_scores;
	}

	public List<Double> boList() {
		// TODO Auto-generated method stub
		return this.bo_scores;
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
			BLAST blast = new BLAST(Oligo.Directory,Oligo.Genome) ;
			HashMap<Integer,String> queries = new HashMap<Integer,String>(); 

			ArrayList< ArrayList<Double>> score_list = new ArrayList< ArrayList<Double> > (this.margin) ;

			// Starting from leftmost position on the span, move to the right and get
			for ( int ii = this.oligo_min ; ii < this.oligo_max ; ii ++) {

				// Extract an oligo for every posible position and put it in the map
				queries.put(ii,SequenceTools.ReverseCompliment(this.getOligo(ii)) );
				//System.err.println(this.getOligo(ii));
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
		catch (Exception ee){
			System.err.println("Error calculating BLAST scores");
			ee.printStackTrace();
			populate_blank_bg();
		}
		
	}

	/**If an error was thrown populating the bg scores, populate the list with all scores set to 0
	 */
	private void populate_blank_bg(){
		try{
			ArrayList< ArrayList<Double>> score_list = new ArrayList< ArrayList<Double> > (this.margin) ;
			
			int count = 0;
			
			for (ArrayList<Double> position : score_list){
				bg_scores.add( count ,0.0);
				count ++;
				bg_sorted.add(count);
			}
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates the Blast Oligo scores for every oligo on the genome given the current valid set
	 * @throws Exception
	 */
	public void calc_bo() throws Exception {

		// Clear all the values in the BO score list and sorted list
		this.bo_scores.clear();
		this.bo_sorted.clear();
		this.scores.clear();

		// Set every possible oligo on the span as the optimized oligo, then calculate the associated BO Value
		for ( int ii = this.oligo_min; ii < this.oligo_max; ii++ ) {
			this.set(ii);
			Double score = mage.Switches.BlastOligo.score(this);
			this.bo_scores.add(score);

			// Add the Index Number for sorting
			this.bo_sorted.add(ii);

			this.scores.add(this.scoreAt(ii));
		}

		// Sort Ascending by BO value
		Collections.sort(this.bo_sorted, new Comparator <Integer>() {

			//Implementing stanard compare fucntion for sorting ascending
			public int compare(final Integer ii, final Integer jj) {
				return (int) ( (bo_scores.get(ii-1))- (bo_scores.get(jj-1)));
			}
		});

		this.greedy_choice = mage.Switches.Oligo.greedyScore(this);
		this.reset();

	}

	/**
	 * Finalize a choice of oligo based on the greedy Score
	 * 
	 */

	public void finalize() throws Exception {

		// Set the oligo to the greedy Choice to recalculate the mistargets
		this.select();

		// Recalculate the score at the final position and update it
		Double score = mage.Switches.BlastOligo.score(this);
		this.bo_scores.set(this.greedy_choice-1, score);
		this.scores.add(this.scoreAt(greedy_choice));

		// Re select this position to update the greedy score;
		this.select();

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
	 * Calculates the BO Score, weighted for the optimized oligo
	 * 
	 * This uses switches.BlastOligoWeight.score 
	 * @throws Exception 
	 * 
	 */
	public void calc_primary_bo() {

		// Calculate and assign the blast oligo score
		try {
			//this.setOligo(this.primary_position);
			this.bo_weighted= mage.Switches.BlastOligo.score(this);
			this.reset();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * 	Calculates the primary score associated with the span.  Use the getPrimaryScoreAsString method to retreive it
	 */
	public void calc_primaryScore() {
		this.primary_position = mage.Switches.PrimaryScore.score(bg_scores,dg_scores) + this.oligo_min;	
	}

	private void calcOptimizedBounds(int start_position){
		this.opt_start	= start_position;
		this.opt_end	= start_position+Oligo.ideal_length-1;
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

	private Double dg(int position){
		if (position >= this.oligo_max && position < this.oligo_min);
		return this.dg_scores.get(position-this.oligo_min);
	}

	public List<Double> dgList() {
		return this.dg_scores;
	}

	public String getBGasString(){	
		StringBuilder sb =  new StringBuilder();

		for (int ii = this.oligo_min, jj=0 ; ii< oligo_max ; ii++, jj++) {
			sb.append(String.format( "%d \t \t %.3f \n",jj, bg(ii)) ); 
		}

		return sb.toString();
	}


	/**
	 * Get A List of Indicies corresponding to sorted BG Positions
	 * @return
	 */
	public List<Integer> getBGSorted() { return this.bg_sorted; }

	public String getBOasString() {

		StringBuilder sb =  new StringBuilder();
		for (int ii = this.oligo_min, jj=0 ; ii< oligo_max ; ii++, jj++) {
			sb.append(String.format( "%d \t \t %.3f \n",ii, this.bo_scores.get(jj) )); //, dg_raw_scores.get(ii-this.oligo_min) )); 
		}

		return sb.toString();

	}

	/**
	 * Get A Lst of Indicies corresponding to sorted BO positions
	 * @return
	 */
	public List<Integer> getBOSorted() {return this.bo_sorted; }

	public String getDGasString(){

		StringBuilder sb =  new StringBuilder();
		for (int ii = this.oligo_min, jj=0 ; ii< oligo_max ; ii++, jj++) {
			sb.append(String.format( "%d \t \t %.3f \n",jj, dg(ii) )); //, dg_raw_scores.get(ii-this.oligo_min) )); 
		}

		return sb.toString();
	}

	/**
	 * Gets a List of Indicices corresponding to sorted DG positions
	 * @return
	 */
	public List<Integer> getDGSorted() { return this.dg_sorted; }

	/**
	 * Returns the last position on the genome
	 * @return
	 */
	public int getGenomeEnd(){	return this.genome_end; }

	/**
	 * Returns the first most position on the genome
	 * @return
	 */
	public int getGenomeStart(){ return this.genome_start; }

	/**
	 * Returns the greedyscore
	 * @return
	 */
	//TODO: Switch to dg_scores?
	private double getGreedyScore() { return bo_scores.get(this.greedy_choice-1); }

	/**
	 * Get the number for possible oligos on the span i.e the length of the margin
	 * @return Integer count of possible oligos
	 */
	public int getMarginLength() {return this.margin;}

	/**
	 * Returns the number of features associated with this oligo
	 * @return	Returns an positive integer representing the number of linked mistarget
	 */
	public int getMistarget(){
		return this.valid_mt.size();
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

	/**
	 * Every oligo has a unique static id number
	 * @return	An integer holding the unique id number of the referenced oligo
	 */
	public int getOligoId(){
		return this.oligo_id;
	}

	/**
	 * Get the last possible starting position on the margin
	 * @return	Integer with starting index
	 */
	public int getOligoMaximum() {return this.oligo_max;}

	/**
	 * Get the first possible starting position on the margin
	 * @return	Integer with starting index
	 */
	public int getOligoMinimum() {return this.oligo_min;}

	/**
	 * @return the oligoType
	 */
	public OligoType getOligoType() {
		return oligoType;
	}

	/**
	 * Returns the optimized sequence as a string of length = Oligo.ideal_length;
	 * @return String of ATCGs
	 */
	public String	getOptimized() { return this.sequence;}

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

	/**
	 * Returns the start position of the oligo that represents the primary position
	 * @return	Integer containing the starting position of the primary oligo.
	 */
	public int getPrimaryPosition() {
		return this.primary_position;
	}

	/**
	 * The primary score is the optimal free energy and blast genome score combined
	 * @return	A formatted string in teh form "Oligo_%d Primary Score ( %.3f , %.3f )"
	 */
	public String getPrimaryScoreAsString(){
		return String.format("Oligo_%d Primary Score ( %.3f , %.3f )",this.oligo_id , this.bg(this.primary_position),this.dg(this.primary_position)); 
	}

	public String getTarget(){
		return this.target;
	}

	/**
	 * Get a list of Indicies corresponding to DG positions below the designated threshold
	 * @return
	 */
	public List<Integer> getValidDG() { return this.valid_dg_positions;}

	public double getWeightedBOScore() {
		return this.bo_weighted;
	}

	/**
	 * Reset the span to make the bounds of the optimized oligo capture the entire span 
	 * @throws Exception	
	 */
	public void reset() throws Exception {
		this.optimized	= sequence;
		this.opt_start 	= this.oligo_min;
		this.opt_end	= this.oligo_max +Oligo.ideal_length;

		this.valid_mt.removeAll(this.valid_mt);
		// For each associated mistarget check if it is valid
		for (Mistarget mt : mt_collection){
			if (mt.isValid(this)) {
				this.valid_mt.add(mt);
			}
		}

		// Reset the oligo score to be N/A for the entire span
		this.currentScore = new OligoScore();

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
	 * Function sets the optmized oligo on the span to be the one defined by the greedy choice
	 * 
	 */
	public void select() {
		try {
			this.set(this.greedy_choice);
			this.currentScore = this.scoreAt(this.greedy_choice);
		} catch (Exception e) {
			e.printStackTrace();
		}
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
	public void set(int new_start_position) throws Exception {
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
	 * Helper Function for setting the optimized bounds
	 * @param start_position
	 * @throws Exception
	 */
	private void setOptimized(int start_position) throws Exception{
		this.optimized = getOligo(start_position);
		calcOptimizedBounds(start_position);
	}

	/**
	 * Return a string containing the nucleotides of the entire span
	 * @return	String of nucleotide A,T,C,G
	 */
	public String getAsString() { return this.sequence; }


	/**
	 * Returns the greedy choice index, remember that this is this position on the margin, indexing begins from 1
	 * @return	Integer between 1 and margin length inclusive
	 */
	public int	getGreedyChoice() { return this.greedy_choice; }

	/**
	 * Returns a list of Oligo Scores
	 * @return	List of OligoScore objects
	 */
	public List<OligoScore> getScores() {return this.scores;}

	/**
	 * Returns the current subsequence selected
	 * @return	A String of ideal oligo length
	 */
	public String getOptimizedAsString() {return this.optimized;}

	/**
	 * Sets the optMagePosition
	 * @param position
	 */
	public void  setOptMagePosition( int position )
	{
		this.optMagePosition = position;
	}

	/**
	 * Returns the position that optMAGE selected
	 * @return
	 */
	public int  getOptMagePosition( )
	{
		//System.err.println("[DEBUG getOptMagePosition] " + this.optMagePosition);
		return this.optMagePosition;
	}
	
	/**
	 * @return the margin
	 */
	public int getMargin() {
		return margin;
	}

	/**
	 * Get a List of all possible oligos covering the valid span
	 * @return
	 * @throws Exception
	 */
	public List<String> getPossibleOligos() throws Exception{

		ArrayList<String> poligos = new ArrayList<String> ();
		for (int ii = this.oligo_min; ii<this.oligo_max; ii++) {
			poligos.add(this.getOligo(ii));
		}
		return poligos;
	}

	/**
	 * This funciton will return the optimal position of shift for the sub-oligo from the span 
	 * with the best score after optimization.
	 * 
	 * @return a 0 indexed number
	 */
	public int getOptimalPosition() {
		return this.getGreedyChoice()-1;
	}

	public static void resetCount() {
		// reset the oligo count
		Oligo.oligo_count = 0;
	}
}
