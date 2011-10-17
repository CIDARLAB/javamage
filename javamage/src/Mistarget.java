import java.util.List;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.template.AbstractSequence;


public class Mistarget implements FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> {

	private String type;
	private String description;
	private SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound> location;
	private Double score;
	private String source;
	private Object userobject;
	private FeatureInterface <AbstractSequence<NucleotideCompound>,NucleotideCompound> parent;
	private List< FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> children;

	public Mistarget(int start, int end, Double score){
		setLocation(start, end);
		setScore(score);
		
		// Ones that done have to be initialized
		this.source= "";
		this.userobject = "";
		this.type = "Homology Site";
		
		// The ones that we don't use...
		this.children = null;
		this.parent = null;
		this.description = "";
		
	}
	
	
	public Double getScore() {
		return this.score;
	}
	
	public void setScore(Double arg0) {
		this.score = arg0;
	}
	
	public void setLocation(int start, int end){
		this.location = new SequenceLocation<AbstractSequence<NucleotideCompound>, NucleotideCompound>(start, end, null);
	}
	
	@Override
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
