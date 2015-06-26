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
import java.util.logging.Level;
import java.util.logging.Logger;
import mage.Tools.Melt;

/**
 *
 * @author mquintin
 */
public class TestMelt {
    public static void main(String[] args) throws IOException, InterruptedException{
        //try {
            callScript();
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
}
