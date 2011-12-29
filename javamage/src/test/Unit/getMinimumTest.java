package test.Unit;

import java.util.ArrayList;



public class getMinimumTest {

	public static void main( String[] args) {
		ArrayList<Double> l1= new ArrayList<Double>();
		ArrayList<Double> l2= new ArrayList<Double>();
		
		l1.add(94.0);
		l1.add(97.0);
		l1.add(94.0);
		l1.add(97.0);
		l1.add(98.0);
		l1.add(14.0);
		l1.add(14.0);
		l1.add(34.0);

		l2.add(0.0);
		l2.add(0.0);
		l2.add(0.0);
		l2.add(0.0);
		l2.add(0.0);
		l2.add(4.0);
		l2.add(1.0);
		l2.add(0.0);
		
		int index = mage.Switches.PrimaryScore.getMinimum(l1,l2);
		System.out.println( index +  " index gives "+ l1.get(index) ) ;
		
	}
	
}
