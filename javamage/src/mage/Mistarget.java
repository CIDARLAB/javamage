package mage;

import tools.BLAST.BlastResult;

public class Mistarget {

	public	int 	a_start;
	public 	int 	a_end;
	public 	int 	b_start;
	public	int 	b_end;
	private	Double 	bitscore;
	private	Double 	evalue;
	public	int 	id_A;
	public	int 	id_B;
	private Oligo	spanA;
	private Oligo	spanB;
	private	String 	sequence;

	private int		a_overlap_start;
	private int		a_overlap_end;
	private int 	b_overlap_start;
	private int		b_overlap_end;

	Mistarget(BlastResult br,int spanA, int spanB) {

		// Record the span IDs and save the oligo reference
		this.id_A   =	spanA; 
		this.id_B   =	spanB;
		this.spanA	=	Oligo.oligo_map.get(this.id_A);		
		this.spanB	=	Oligo.oligo_map.get(this.id_B);

		// Get the start and end values in the query
		this.a_end  =	br.sEnd;
		this.a_start= 	br.sStart;
		this.b_end 	= 	br.qEnd;
		this.b_start= 	br.qStart;

		// Assign bitscore and evalue of the mistarget
		this.bitscore=	br.bitscore;
		this.evalue  =	br.evalue;

		// Set overlaps by default to lengths
		this.a_overlap_start= 0;
		this.b_overlap_start= 0;
		this.a_overlap_end	= sequence.length();
		this.b_overlap_end	= sequence.length();

		// Register this mistarget with those spans
		this.spanA.addMistarget(this);
		this.spanB.addMistarget(this);

		System.err.println("Mistarget Between Oligo "+id_A+" and Oligo "+id_B+" : "+this.sequence);

	}


//	public Double getWeightedScore() {
//		//Switches.WeightOverlap.Score();
//	}
	
	/**
	 * Calculates the number of overlapping basepairs for the given spans and weights accordingly
	 * 
	 * @return	Returns integer count of overlapping basepairs
	 */
	public int overlap(){

		int count =0;
		// Bound checking for the entire sequence. O(n) complexity and but with n < 1000 should be trivial w/primatives
		for (int ii=0; ii<sequence.length();ii++){
			if ( 	(ii <= a_overlap_end && ii >= a_overlap_start ) 
					&&	(ii <= b_overlap_end && ii >= b_overlap_start )	)
			{count++;}
		}

		return count;
	}



	/**
	 * 
	 * isValid will checks the position of the mistarget on the oligo and returns true if any part of the mistarget is found on the span
	 * 
	 * @param oligo				The span associated with the mistarget
	 * @param start_position	The start position of the oligo on the span
	 * @param end_position		The end position of the oligo on the span
	 * @return					Returns true of the mistarget is still on the span, false if the mistarget is not on the span
	 * @throws Exception		In the event that you passed an oligo that is not associated with this mistarget, then the oligo is removed
	 */
	public boolean isValid(Oligo oligo, int start_position, int end_position) throws Exception {

		int mstart	= this.a_start; 
		int mend	= this.a_end;

		int overlap_start = -1;
		int overlap_end	  = -1;

		boolean valid 	= true;
		boolean isA		= true;
		// if the oligo in reference is not A, we assume B
		if (oligo.getOligoId() == this.id_B )
		{
			isA		= false;
			mstart	= this.b_start;
			mend	= this.b_end;
		}
		else if ( oligo.getOligoId() != this.id_A) {
			// This situation should never occur, if it does that means that there is an error in the way the mistargets are refereneced amongst their relative oligos
			throw new Exception("[Mistarget] Fatal Error - Mismarget Referenced Incorrectly");
		}

		// Does the mistarget start after the oligo ends?
		if (mstart > end_position)	{ valid = false;}
		// Does the mistarget end before the oligo starts?
		if (mend < mstart)			{ valid = false;}

		// Check if we have 100% overlap
		if (mend <= end_position && mstart >= start_position )	{ overlap_start=0; overlap_end=this.sequence.length(); }

		// Check if we have mistarget hanging of the end
		if (mend > end_position && mstart >= start_position)	{ overlap_start=0; overlap_end=(end_position-mstart); }

		// Check if we have the mistartget hanging of the start
		if (mend <= start_position && mstart < start_position)	{ overlap_start= (mend - start_position +1) ; overlap_end=this.sequence.length(); }

		// Assign it to the correct variable
		if (isA){
			a_overlap_start = overlap_start;
			a_overlap_end	= overlap_end;
		}
		else {
			b_overlap_start	= overlap_start;
			b_overlap_end	= overlap_end;
		}

		// if the above criteria is not met, return true
		return valid;
	}


}