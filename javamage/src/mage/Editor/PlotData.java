package mage.Editor;

import mage.Core.Oligo;

/**
 * This is an object for storing the plot values of an Oligo
 * Intended to passed onto mage editor.
 * 
 * @author mockingbird
 */
public class PlotData {

	private Double[] FreeEnergy;
	private Double[] BlastGenome;
	private Double[] BlastOligo;
	private int optMagePosition;
	private int merlinPosition;
	
	/**
	 * Create Mage-Editor Plot Object
	 * This object holds data for plotting in mage editor
	 * 
	 * @param ol Given an oligo, this object will return a plot object
	 */
	public PlotData( Oligo ol ){
		
		// Create Double Arrays
		this.FreeEnergy = ol.dgList().toArray(new Double[ol.dgList().size()]);
		this.BlastOligo = ol.boList().toArray(new Double[ol.boList().size()]);
		this.BlastGenome = ol.bgList().toArray( new Double[ol.bgList().size()] ) ;
		
		// Extract optMageposition and merlin position
		this.optMagePosition = ol.getOptMagePosition();
		this.merlinPosition = ol.getOptimalPosition();
	}
	
	/**
	 * Get the free Energy score as a string
	 * @return free energy score string
	 */
	public String getFreeEnergy()
	{
		return getArrayAsString(this.FreeEnergy);
	}
	
	/**
	 * Get the Blast Oligo score as a string
	 * @return	BO Score as string
	 */
	public String getBlastOligo(){
		return getArrayAsString(this.BlastOligo);
	}
	
	/**
	 * Get the bo score as a string
	 * @return bo score
	 */
	public String getBlastGenome(){
		return getArrayAsString(this.BlastGenome);
	}
	
	/**
	 * Get the shift position of the corresponding optimized oligo as calculated by optmage
	 * @return	A number in the form of a string
	 */
	public String getOptMage(){
		return Integer.toString(this.optMagePosition);
	}
	
	/** Get the position that was selected by merlin
	 * 
	 * @return string with integer position of merlin selected
	 */
	public String getMerlin(){
		return Integer.toString(this.merlinPosition);
	}
	
	/**
	 * Helper function, given an array of doubles, this will return a string
	 * in which each element of the input array is on a new line
	 * @param array  Array of Double values to be converted to String
	 * @return Single string with each element separated by newlines
	 */
	private String getArrayAsString(Double [] array)
	{
		StringBuilder sb = new StringBuilder();
		for ( Double dd : array )
		{
			sb.append(dd.toString()+"\n");
		}
		return sb.toString();
	}
}
