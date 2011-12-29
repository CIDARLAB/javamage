package test.Unit;

import java.io.IOException;

import mage.Core.Oligo;
import test.Constants;

public class FASTA {

public static void main(String[] args) throws IOException {
		
		System.out.println(mage.Tools.FASTA.readFFN(Constants.blastdirectory,Oligo.Genome).length() );
	}
	
}
