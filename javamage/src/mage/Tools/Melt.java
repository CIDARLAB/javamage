package mage.Tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A wrapper for merlinmelt.py on a list of oligos and collecting 
 * the results. The script uses on BioPython to calculate Tm_NN
 * @author Michael Quintin, adapted from code by Samir Ahmed
 *
 */
public abstract class Melt {
    
    protected static String result = "";
    protected static String script = "merlinmelt.py";
    //private static final String separator = java.io.File.separator;
    //protected static String scriptpath = "\\.." + separator;

    /**
     * 
     * @param args
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    public static String execute(String args) throws IOException, InterruptedException {   
        //find the script
        ClassLoader loader = Melt.class.getClassLoader();
        URL url = loader.getResource(script);
        System.out.println(url.toExternalForm());
        System.out.println(url.toString());
        String scriptPath = url.getFile();
        if (scriptPath.startsWith("/") | scriptPath.startsWith("\\")){
            scriptPath = scriptPath.substring(1);
        }
        
        args = Constants.python + " " +  scriptPath +  " " + args;
        String[] argarr = args.split(" ");
        ProcessBuilder pb = new ProcessBuilder(argarr); //yes, I know we merged
            //the args into one string earlier. Consider this whitespace handling
        Process p = pb.start();
        p.waitFor();

        BufferedReader reader = new BufferedReader(
                new InputStreamReader(p.getInputStream()));
        StringBuilder builder = new StringBuilder();
        String line = null;
        while ( (line = reader.readLine()) != null) {
           builder.append(line);
           builder.append('\t');
        }
        result = builder.toString();
        return result;
    }
    
    //convert the arguments into a more friendly format to be passed to the script
    public static String execute(String[] args) throws IOException, InterruptedException{
        String s = "";
        for (String arg : args){
            s = s.concat(arg).concat(" ");
        }
        return execute(s);
    }
    public static String execute(List<String> args) throws IOException, InterruptedException{
        String[] s = args.toArray(new String[args.size()]);
        return execute(s);
    }
    
    
}
