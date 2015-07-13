/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test.Unit;

/**
 *
 * @author mquintin
 */
public class TestReverseCompliment {
    
    public static void main(String[] args){
        String s = "CAT";
        String expect = "ATG";
        String res = mage.Tools.SequenceTools.ReverseCompliment(s);
        System.out.println(expect.equals(res));
        
        System.out.println(mage.Tools.SequenceTools.ReverseCompliment(s));
    }
    
}
