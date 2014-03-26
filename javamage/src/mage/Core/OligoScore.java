package mage.Core;


/**
 * 
 * Class for keeping track of the score of the oligo
 * 
 * @author Samir Ahmed
 *
 */
public class OligoScore {

	private	Double	BO; //Blast Oligo
	private	Double	DG; //Free Energy
	private	Double	BG; //Blast Genome
	
	
	/**
	 * Construct a instance in which all the components of the score are entirely NaN
	 * Use the getters and setters to define each component individually
	 */
	public OligoScore() {
		this.BG = Double.NaN;
		this.DG = Double.NaN;
		this.BO = Double.NaN;
	}
	
	/**
	 * Construct a new instance of a Global score
	 * @param BG	The Blast Genome Component
	 * @param DG	The Free Energy Component
	 * @param BO	The Blast Oligo Component
	 */
	public OligoScore(double BG, double DG, double BO){
		this.BO = BO;
		this.BG = BG;
		this.DG = DG;
	}
	
	/**
	 * Returns true if the new score is better than the old score
	 * 
	 * @param newScore
	 * @return
	 */
	public boolean isBetterThan(OligoScore bestScore) {
		return ( mage.Switches.Oligo.compare(bestScore,this) > 0 ) ? true : false ;
	}
	
	
	/**
	 * Get Blast Oligo Component
	 * @return	Total Blast Oligo Score
	 */
	public Double BlastOligo() {
		return this.BO;
	}
	
	/**
	 * Get Blast Genome Compoenent
	 * @return	Total Blast Genome Score
	 */
	public Double BlastGenome() {
		return this.BG;
	}
	
	/**
	 * Get Free Energy Component
	 * @return	Total Free Energy Score
	 */
	public Double FreeEnergy() {
		return this.DG;
	}
	
	/**
	 * Set the Blast Oligo Score Component
	 * @param newBO	New Value
	 */
	public void BlastOligo(double newBO){
		this.BO = newBO;
	}
	
	/**
	 * Set the Blast Oligo Score Component
	 * @param newBG	New BG Value
	 */
	public void BlastGenome(double newBG){
		this.BG = newBG;
	}
	
	
	/**
	 * Set the Free Energy Component
	 * @param newDG	New DG value
	 */
	public void FreeEnergy(double newDG) {
		this.DG = newDG;
	}
	
	/**
	 * Returns a formatted string
	 */
	public String toString(){
		return String.format("{%7.3f,%7.3f,%7.3f}", this.DG,this.BG, this.BO);
	}

	/**
	 * Will add the score values of the 
	 * 
	 * @param currentScore
	 */
	public void add(OligoScore currentScore) {
		mage.Switches.Oligo.add(this,currentScore);
	}
	
}
