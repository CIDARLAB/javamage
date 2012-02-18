package mage.Core;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import utils.TextFile;

/**
 * 
 * @author Samir Ahmed
 *
 */
public class Target {

	final public	String	gene_name;
	final public 	boolean	sense;
	final public	int	 	replichore;
	final public	int		left_position;
	final public	int		right_position;
	final public	String	type;
	final public	String	sequence;

	/**
	 * Constructs a target from a the given line of an INPUTtarg.txt file
	 * @param line
	 */
	public Target(String line){
		String [] properties = line.split("\t");
		this.gene_name 		= properties[0];
		this.sense	   		= properties[1].startsWith("+") ? true : false;
		this.replichore		= Integer.parseInt(properties[2]);
		this.left_position	= Integer.parseInt(properties[3]);
		this.right_position	= Integer.parseInt(properties[4]);
		this.type			= properties[5];
		
		if ( properties.length > 6 ){ this.sequence	= properties[6]; }
		else {	this.sequence= ""; }
	}

	/**
	 * Creates a list a targets from a given file
	 * @param 	filepath	File path of the INPUt target file
	 * @return	List of desired targets
	 * @throws 	IOException
	 */
	public static List<Target> loadTarget(String filepath) throws IOException{
		String inputTargetFile 	= TextFile.read(filepath);
		String[] lines 		= inputTargetFile.split("\n");

		ArrayList<Target> list	= new ArrayList<Target>(lines.length-1);
		for (int ii=1; ii<lines.length; ii++) {
			if (lines[ii]!=null && !lines[ii].isEmpty()) {
				Target tt =  new Target(lines[ii]);
				list.add(tt);
			}
		}

		return list;
	}

	/**
	 * Creates a list of targets from a given file
	 *  
	 * @param directory	File directory
	 * @param filename	File name
	 * @return			A list of Targets
	 * @throws IOException
	 */
	public static List<Target> loadTarget(String directory,String filename) throws IOException{
			return loadTarget(directory+filename); 
	}

}
