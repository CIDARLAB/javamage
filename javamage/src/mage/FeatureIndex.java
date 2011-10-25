package mage;


/**
 * Simple Index Holder
 * This allows us to keep track on every oligo feature by oligo parent number, and child index number
 * 
 * @author Samir Ahmed
 *
 *	e.g  New FeatureIndex(10,3);
 *	this means that this is the feature (homology site) number 3 on Oligo number 10
 *
 */
public class FeatureIndex {

	public int parent_oligo;
	public int child_index;

	/**
	 * Constructor for a FeatureIndex
	 * 
	 * @param Oligo The Index Number of the parent Oligo
	 * @param Index The Index of the possessing feature
	 */
	FeatureIndex(int Oligo, int Index){
		this.parent_oligo = Oligo;
		this.child_index = Index;
	}
}
