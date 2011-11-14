package mage.Switches;
//
//import java.util.ArrayList;
//
//import mage.Oligo;
//
//public class FullOligo {
//
//	public static int method = 1;
//
//	public static ArrayList<Double> score(Oligo ol) throws Exception {
//
//		ArrayList<Double> bo_scores = new ArrayList<Double>(ol.getMarginLength());
//
//		switch (mage.Switches.FullOligo.method) {
//
//		// Just sum all the mistarget raw scores.
//		case 1:	
//
//			for ( int ii = ol.getOligoMinimum(); ii < ol.getOligoMaximum(); ii++ ) {
//				ol.setOligo(ii);
//				Double score = mage.Switches.BlastOligo.score(ol);
////				for (Mistarget mt : ol.valid_mt) {
////					score += mt.score();
////				}
//				bo_scores.add(score);
//			}
//
//			//		for (int ii = ol.getOligoMinimum() )
//			//		for(Mistarget mt : ol.valid_mt) { score += mt.getRawScore(); } break;
//		default : break;
//		}
//
//		return bo_scores;
//
//	}
//}