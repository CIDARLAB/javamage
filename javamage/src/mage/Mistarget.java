package mage;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.template.AbstractSequence;

import tools.BLAST.BlastResult;


public class Mistarget { //implements FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> {

	private String type;
	private FeatureIndex index;
	private String description;
	private SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound> location;
	private Double score;
	private String source;
	private Object userobject;
	private List< FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> children;
	
//	private Double = raw_score;
	
	private int parent;
	private String sequence;
	private ArrayList<Double> raw_list;
	private ArrayList<Double> weighted_list;
	private HashSet<Integer> connected_oligos;
	private HashMap<Integer,FeatureIndex> location_map;
	private HashMap<Integer,Double> raw_map;
	private HashMap<Integer,Double> weighted_map;
	private Integer start;
	private Integer end; 
	public static HashMap<FeatureIndex,Integer> mistarget_index;

	public Mistarget(BlastResult br, int subkectOligo, int parentOligo) {
		//this(br,subjectOligo,parentOligo,false);
	}
	
	public Mistarget(BlastResult br, int subjectOligo, int parentOligo, boolean IsNew){
		
		// Given the input arguments, we can distinguish between parent and linked Oligos
		// After distinguishing the Oligos, we create an identical Mistarget on for the linked Oligo
		
		// Get the parent Oligo value and store the sequence value
		this.parent = parentOligo;
		this.sequence = br.match_sequence;
		
		// Determine the connected Oligo
		int otherID;
		if (parentOligo ==  subjectOligo) {
			this.start = br.sStart;
			this.end = br.sEnd;
			otherID = br.oligoID;
		} 
		else {
			this.start = br.qStart;
			this.end = br.qEnd;
			otherID = subjectOligo;
		}
		
		// For the first mistarget of this sequence 
		if (IsNew) {
			//Mistarget.addOligo(OligoID, br, subjectOligo, otherID);
		}
		
		// Score the Oligo
		
		// Weight the score
		
		// 
		
		// Create a linked mistarget on the other Oligo;
		
		
		// Add the connection
		
		// Initialize some of the arrays
		
		
	
		//addConnection();
				
	}
	
	
//	private addScore();
//	
//	public void addOligo(BlastResult br,int subjectID){
//		
//		if ( !(connected_oligos.contains(br.oligoID) && connected_oligos.contains(subjectID) ) )
//		{
//			if (!connected_oligos.contains(br.oligoID)) {
//				connected_oligos.add(br.oligoID);
//				location_map.put(subjectID,new FeatureIndex(br.qStart,br.qEnd) );
//				raw_list.
//			}
//			
//			if (!connected_oligos.contains(subjectID) ) {
//				connected_oligos.add(subjectID);
//				location_map.put(subjectID,new FeatureIndex(br.sStart,br.sEnd) );
//			}
//			
//			Switches.SumMistargetScores(this.raw_scores);
//			
//		}
//		else{throw new Exception("[Mistarget] Mistarget already Exists: "+subjectID + " "+br.oligoID);}
//	}
//	
//		this.raw_score += Switches.BlastScore(br.bitscore , br.evalue);
//		br.oligoID;
//		
//		br.match_sequence;
//		br.qEnd;
//		br.sEnd;
//		br.qStart;
//		br.sStart;		
//	}
//	
	/**
	 * 
	 * Constructor takes in location and scoring information
	 * 
	 * @param index Contains information about the parent oligo and feature number with relation to the feature oligo
	 * @param start	Start location of the mistarget
	 * @param end	End location of the mistarget
	 * @param score	Score of the mistarget as calulated by Switches.MistargetScore
	 */
	
	/*

	public Mistarget(FeatureIndex index ,int start, int end, Double score){
		setLocation(start, end);
		setScore(score);

		// Ones that done have to be initialized
		this.index = index;
		this.source= "";
		this.userobject = "";
		this.type = "Homology Site";

		// The ones that we don't use...
		this.children = null;
		this.parent = null;
		this.description = "";

	}

	//A Set of Constructors that all for a variety of constructor calls
	
	public Mistarget(FeatureIndex index, int start, int end, Double bitscore, Double evalue) {
		this(index, start, end, Switches.BlastScore(bitscore, evalue));
	}

	public Mistarget(int parent, int feature_no, int start, int end, Double score){
		this (new FeatureIndex (parent, feature_no), start, end, score);
	}

	public Mistarget(int parent, int feature_no, int start, int end, Double bitscore, Double evalue){
		this( new FeatureIndex(parent, feature_no), start, end, Switches.BlastScore(bitscore, evalue));
	}

	// Returns the score of a given oligo
	public Double getScore() {
		return this.score;
	}

	// Set the score of a given oligo
	public void setScore(Double arg0) {
		this.score = arg0;
	}

	// Set the score of a given oligo with bitscore and evalue
	public void setScore(Double bitscore, Double evalue) {
		this.score = Switches.BlastScore(bitscore, evalue);
	}

	// Set the Location with start and ent point
	public void setLocation(int start, int end){
		this.location = new SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound>(start, end, null);
	}
	
	// Return the Feature Index Location about the parent
	public FeatureIndex getFeatureIndex (){
		return this.index;
	}
	
	// Get just the Parent Oligo Information
	public int getParentOligo() {
		return this.index.parent_oligo;
	}
	
	// Get just the Child Feature Index number Information
	public int getChildIndex(){
		return this.index.parent_oligo;
	}

	@Override
	// Set the Location with a Sequence Location Object
	public void setLocation( SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound> arg0) {
		this.location =  arg0;
	}

	// Everything else is for implementing the interface

	@Override
	public SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound> getLocations() {
		return this.location;
	}


	@Override
	public String getSource() {
		return this.source;
	}

	@Override
	public String getType() {
		return this.type;
	}

	@Override
	public Object getUserObject() {
		return this.userobject;
	}

	@Override
	public void setSource(String arg0) {
		this.source = arg0;
	}

	@Override
	public void setUserObject(Object arg0) {
		this.userobject = arg0;
	}


	@Override
	public List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> getChildrenFeatures() {
		return this.children; 
	}

	@Override
	public String getDescription() { 
		return this.description; 
	}

	@Override
	public FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> getParentFeature() {
		return this.parent;
	}

	@Override
	public String getShortDescription() {
		return this.description;
	}

	@Override
	public void setDescription(String arg0) {}

	@Override
	public void setShortDescription(String arg0) {}

	@Override
	public void setType(String arg0) {}

	@Override
	public void setChildrenFeatures(
			List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> arg0) {}

	@Override
	public void setParentFeature(
			FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> arg0) {}
*/
}
