package test.Unit;

import mage.Core.Target;
import test.Constants;

import java.util.List;



public class TargetTesting {

	
	// Test function for loading Target and each property
	public static void main(String[] args) throws Exception{
		
		List<Target> list = Target.loadTarget(Constants.ioDirectory,"INPUTtarg.txt");
		
		for (Target tt : list){
			
			System.out.print(tt.gene_name);
			System.out.print("\t"+tt.sense);
			System.out.print("\t"+tt.replichore);
			System.out.print("\t"+tt.left_position);
			System.out.print("\t"+tt.right_position);
			System.out.print("\t"+tt.type);
			System.out.println("\t"+tt.sequence);
		
		}
	}
}
