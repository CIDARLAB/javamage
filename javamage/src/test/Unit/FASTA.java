package test.Unit;

import java.io.IOException;

import mage.Core.Oligo;
import test.Constants;
import utils.TextFile;

public class FASTA {

public static void main(String[] args) throws IOException {
		
		String genome =(mage.Tools.FASTA.readFFN(Constants.blastdirectory,Oligo.Genome));
				mage.Tools.FASTA.writeFFN(Constants.blastdirectory,"testwrite",genome.substring(0,1000) );
		System.out.println(TextFile.read(Constants.blastdirectory+"testwrite"));
	}
	
}
