package mage;

import java.util.List;


import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.template.AbstractSequence;


public class Mistarget implements FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> {

	private String type;
	private FeatureIndex index;
	private String description;
	private SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound> location;
	private Double score;
	private String source;
	private Object userobject;
	private FeatureInterface <AbstractSequence<NucleotideCompound>,NucleotideCompound> parent;
	private List< FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> children;

	/**
	 * 
	 * Constructor takes in location and scoring information
	 * 
	 * @param index Contains information about the parent oligo and feature number with relation to the feature oligo
	 * @param start	Start location of the mistarget
	 * @param end	End location of the mistarget
	 * @param score	Score of the mistarget as calulated by Switches.MistargetScore
	 */
	
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

	/* A Set of Constructors that all for a variety of constructor calls*/
	
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

	/* Everything else is for implementing the interface*/

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

	/* UNUSED METHODS */

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

}
