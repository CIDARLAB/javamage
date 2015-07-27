/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test.Unit;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import mage.Tools.Pcr.Melt;

/**
 *
 * @author mquintin
 */
public class TestMelt {
    public static void main(String[] args) throws IOException, InterruptedException{
        //try {
            callScript();
            //parse();
        //} catch (Exception ex) {
        //    System.out.print(ex.getMessage());
        //}
    }

    
    private static void callScript() throws IOException, InterruptedException{
        ArrayList<String> list = new ArrayList<>();
        list.add("atcgatcggctaggtaacagattaatctctcggagctgatacgac");
        list.add("TTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACGGGATTATTATTGGGCTCGAATCTACCGTCGATATTGCTGAGTCCACCC");		
        
        String res = Melt.execute(list);
        System.out.println(res);
    }
    
    private static void parse(){
        String s = " 24    -9 2.8 ";
        Double[] expect = new Double[3];
        expect[0] = Double.valueOf(24);
        expect[1] = Double.valueOf(-9);
        expect[2] = Double.valueOf(2.8);
        Boolean res = Arrays.equals(expect,Melt.parseResults(s));
        System.out.println("Melt.parseResults correct: " + res);
        if (!res){
            System.out.println("Expected " + expect);
            System.out.println("Got " + Melt.parseResults(s));
        }
    }
}
