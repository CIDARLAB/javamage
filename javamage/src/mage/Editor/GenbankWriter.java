package mage.Editor;

import utils.TextFile;
import mage.Core.Oligo;


/**
 * 
 * GenbankWriter is a wrapper for creating a genbank file from Merlin Results.
 * The features are intended to be read by MageEditor
 * 
 * @author SamirAhmed
 *
 */
public class GenbankWriter {

	final public static String LOCUS = "LOCUS";
	final public static String ORIGIN = "ORIGIN";
	final public static String FEATUES = "FEATURES";
	final public static String KEYWORDS = "KEYWORDS";
	final public static String ACCESSION = "ACCESSION";
	final public static String VERSION = "VERSION";
	final public static String ENDFILE = "//\n";
	final public static String BP = "bp";
	final public static String LINEAR = "linear";
	final public static String LABEL = "/label=";
	final public static String BO = "/BO=";
	final public static String BG = "/BG=";
	final public static String DG = "/DG=";
	final public static String UNKNOWN = "unknown";
	final public static String Date = "01-JAN-2000";
	final public static String LOCATION_TAG = "Location/Qualifiers";
	final public static String MERLIN = "merlin";
	final public static String MERLIN_SELECTION = "Merlin";
	final public static String OPTMAGE = "optmage";
	final public static String OPTMAGE_SELECTION = "Optmage";
	private static final int CHUNK_COUNT = 6;
	private static final String MUTATION = "mutation";
	private static final String NEW = "New Base Pairs";
	
	private String name ;  
	private String length ;
	private String merlinBO;
	private String merlinBG;
	private String merlinDG;
	private String optMageBO;
	private String optMageBG;
	private String optMageDG;
	private String optMageStart;
	private String optMageStop;
	private String merlinStart;
	private String merlinStop;
	private StringBuilder sb;

	/**
	 * Given an oligo, a genbank file can made with specific parameters for passing back to mage-editor
	 * 
	 * @param ol Oligo for which to generate a genbank file
	 */
	public GenbankWriter(Oligo ol)
	{
		// Setup the data members
		//System.err.println("[DEBUG GENBANK] setup the data members");
		this.name = ol.name;
		this.length = Integer.toString(ol.span);
		
		// Generate the first 5 rows
		//System.err.println("[DEBUG GENBANK] generate the first five rows");
		this.sb = new StringBuilder();
	    this.sb.append(String.format("%-12s%s %s %s %s %s\n",GenbankWriter.LOCUS,this.name,this.length,GenbankWriter.BP,GenbankWriter.LINEAR,GenbankWriter.Date));
	    this.sb.append(String.format("%-12s%s\n",GenbankWriter.ACCESSION,GenbankWriter.UNKNOWN));
	    this.sb.append(String.format("%-12s%s\n",GenbankWriter.VERSION,GenbankWriter.UNKNOWN));
	    this.sb.append(String.format("%-12s\n", GenbankWriter.KEYWORDS));
	    this.sb.append(String.format("%-12s%9s%s\n", GenbankWriter.FEATUES,"",GenbankWriter.LOCATION_TAG));
	    
	    // Now add merlin Tags, optmage Tag and mutationTags
		//System.err.println("[DEBUG GENBANK] now add merlin tags, optmage tag and mutation tags");
		//System.err.println("[DEBUG GENBANK] getOptMageInfo");
	    this.sb.append( getOptMageInfo(ol) );
	    //System.err.println("[DEBUG GENBANK] getMerlinInfo");
	    this.sb.append( getMerlinInfo(ol) );
	    
	    // Add a feature to denote any new basepairs in hte sequence.
		//System.err.println("[DEBUG GENBANK] add a feature to denote any new basepairs in the sequence");
	    this.sb.append(getMutationInfo(ol));
	    
	    // Final Add the sequence
		//System.err.println("[DEBUG GENBANK] add the sequence");
	    this.sb.append(String.format("%-12s\n",GenbankWriter.ORIGIN));
	    this.sb.append(genbankSequenceFormat(ol.sequence));
	    
	    // Append the ENDFILE Mark
		//System.err.println("[DEBUG GENBANK] append the endfile mark");
	    this.sb.append(GenbankWriter.ENDFILE);
	}
	
	/**
	 * Given a sequence, creates genbank formated sequence strings
	 * @param sequence
	 * @return	String with sequence in genbank form
	 */
	private String genbankSequenceFormat(String sequence)
	{
		// Get chunks of size 10
		String [] chunks = TextFile.splitEqually(sequence, 10);
		
		StringBuilder sb = new StringBuilder();
		
		for ( int start = 0; start < chunks.length ; start += GenbankWriter.CHUNK_COUNT, sb.append("\n") )
		{
			// Put in the line count
			sb.append(  String.format("%10s", Integer.toString(start*10+1)) );
			
			// For upto CHUNK_COUNT
			for ( int jj = start; jj < Math.min(chunks.length, start+GenbankWriter.CHUNK_COUNT) ; jj++ )
			{
				sb.append(String.format(" %-10s", chunks[jj] ));
			}
		}
		
		return sb.toString();
		
	}
	
	
	private String getOptMageInfo(Oligo ol)
	{
		this.optMageStart = Integer.toString( ol.getOptMagePosition()+1 );
		this.optMageStop = Integer.toString( ol.getOptMagePosition()+Oligo.ideal_length );
		
		//TODO: Confirm this works in the case where OptMage didn't work, so OptMagePosition=-1.
		//so, replace it with 0
		//TODO: Fix.Confirm this doesn't break the scoring process. It should be able to tell when this happened and ignore scores
		//System.err.println("[DEBUG OPTMAGE] get optMageBG");
		this.optMageBG = Double.toString(ol.bgList().get(Math.max(0,ol.getOptMagePosition())));
		//System.err.println("[DEBUG OPTMAGE] get optMageDG");
		this.optMageDG = Double.toString(ol.dgList().get(Math.max(0,ol.getOptMagePosition())));
		//System.err.println("[DEBUG OPTMAGE] get optMageBO");
		this.optMageBO = Double.toString(ol.boList().get(Math.max(0,ol.getOptMagePosition())));
		
		return String.format("%5s%-17s%s..%s\n"+
					"%21s%s\"%s\"\n" +
					"%21s%s%s\n" +
					"%21s%s%s\n" +
					"%21s%s%s\n",
				"", GenbankWriter.OPTMAGE, this.optMageStart, this.optMageStop, 
				"",	GenbankWriter.LABEL,GenbankWriter.OPTMAGE_SELECTION,
				"", GenbankWriter.DG, this.optMageDG,
				"", GenbankWriter.BG, this.optMageBG,
				"", GenbankWriter.BO, this.optMageBO);
		
	}
	
	private String getMerlinInfo(Oligo ol)
	{
		this.merlinStart = Integer.toString( ol.getOptimalPosition()+1 );
		this.merlinStop = Integer.toString(  ol.getOptimalPosition()+Oligo.ideal_length );
		this.merlinBG = Double.toString(ol.bgList().get(ol.getOptimalPosition()));
		this.merlinDG = Double.toString(ol.dgList().get(ol.getOptimalPosition()));
		this.merlinBO = Double.toString(ol.boList().get(ol.getOptimalPosition()));
	
		return String.format("%5s%-17s%s..%s\n"+
					"%21s%s\"%s\"\n" +
					"%21s%s%s\n" +
					"%21s%s%s\n" +
					"%21s%s%s\n",
				"", GenbankWriter.MERLIN, this.merlinStart, this.merlinStop, 
				"",	GenbankWriter.LABEL,GenbankWriter.MERLIN_SELECTION,
				"", GenbankWriter.DG, this.merlinDG,
				"", GenbankWriter.BG, this.merlinBG,
				"", GenbankWriter.BO, this.merlinBO);
	
	}
	
	private String getMutationInfo(Oligo ol)
	{
		if (ol.target_length == 0)
		{
			return "";
		}
		else
		{
			String start =Integer.toString(ol.target_position+1);
			String stop  =Integer.toString(ol.target_position+ol.target_length);
			return String.format("%5s%-17s%s..%s\n"+
								 "%21s%s\"%s\"\n",
								 "",GenbankWriter.MUTATION, start,stop,
								 "",GenbankWriter.LABEL, GenbankWriter.NEW );
		}
		
	}
	
	
	/**
	 * Get the genbank file as a string
	 * 
	 */
	public String toString() {
		return this.sb.toString();
	}
	
	
}
