package mage.Tools.Pcr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import mage.Core.Oligo;
import mage.Core.Primer;
import mage.Tools.Constants;

/**
 * A wrapper for merlinmelt.py on a list of oligos and collecting the results.
 * The script uses on BioPython to calculate Tm_NN
 *
 * @author Michael Quintin, adapted from code by Samir Ahmed
 *
 */
public abstract class Melt {
    //protected static int calls = 0;

    protected static String result = "";
    protected static String script = "merlinmelt.py";
    //private static final String separator = java.io.File.separator;
    //protected static String scriptpath = "\\.." + separator;

    //command-line parameters to pass to the script
    protected static String nn_table;
    protected static String tmm_table;
    protected static String imm_table;
    protected static String de_table;
    protected static Double dnac1 = 25.0; //nM
    protected static Double dnac2 = 25.0;
    protected static Double na = 50.0; //ions are mM
    protected static Double k = 0.0;
    protected static Double tris = 0.0;
    protected static Double mg = 3.0; //TODO: reset to 0
    protected static Double dntps = 0.8; //mM //TODO: reset to 0
    protected static Double saltcorr;

    public static String directory;

    /**
     *
     * @param args
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    public static String execute(String args) throws IOException, InterruptedException {
        String scriptPath = Oligo.Directory + script;
        //confirm it exists
        File f = new File(scriptPath);
        if (!f.isFile()) {
            //depending on the build setup, it might have gone here
            scriptPath = Oligo.Directory + "scripts" + System.getProperty("file.separator") + script;
        }

        args = Constants.python + " " + scriptPath + " " + args + getCmdArgs();
        String[] argarr = args.split(" ");
        ProcessBuilder pb = new ProcessBuilder(argarr); //yes, I know we merged
        //the args into one string earlier. Consider this whitespace handling
        Process p = pb.start();
        p.waitFor();

        BufferedReader reader = new BufferedReader(
                new InputStreamReader(p.getInputStream()));
        StringBuilder builder = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null) {
            builder.append(line);
            builder.append('\t');
        }
        result = builder.toString();
        reader.close();
        return result;
    }

    //convert the arguments into a more friendly format to be passed to the script
    public static String execute(String[] args) throws IOException, InterruptedException {
        String s = "";
        for (String arg : args) {
            s = s.concat(arg).concat(" ");
        }
        return execute(s);
    }

    public static String execute(List<String> args) throws IOException, InterruptedException {
        String[] s = args.toArray(new String[args.size()]);
        return execute(s);
    }

    public static Double[] parseResults(String res) {
        res = res.trim().replaceAll("\\s+", " ");
        String[] sarr = res.split(" ");
        Double[] darr = new Double[sarr.length];
        for (int i = 0; i < darr.length; i++) {
            darr[i] = Double.valueOf(sarr[i]);
        }
        return darr;
    }

    public static Double[] getMT(String[] primers) throws IOException, InterruptedException {
        String res = execute(primers);
        return parseResults(res);
    }

    public static Double[] getMT(List<String> primers) throws IOException, InterruptedException {
        String res = execute(primers);
        return parseResults(res);
    }

    /**
     * for each primer, set its mt property.
     *
     * @param primers
     * @return
     */
    public static void setMTs(ArrayList<Primer> primers) {
        try {
            ArrayList<Primer> hasMT = new ArrayList();
            ArrayList<Primer> needsMT = new ArrayList();
            for (Primer p : primers) {
                if (p.mt == null) {
                    needsMT.add(p);
                } else {
                    hasMT.add(p);
                }
            }

            String[] seqs = new String[needsMT.size()];
            for (int i = 0; i < needsMT.size(); i++) {
                seqs[i] = needsMT.get(i).seq;
            }
            Double[] mts = getMT(Arrays.asList(seqs));
            for (int i = 0; i < needsMT.size(); i++) {
                needsMT.get(i).setMt(mts[i]);
            }

        } catch (IOException | InterruptedException ex) {
            Logger.getLogger(Melt.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void setMTs(Primer[] primers) {
        ArrayList<Primer> arr = new ArrayList(Arrays.asList(primers));
        setMTs(arr);
    }

    public static void setCmdArg(String name, String value) {
        switch (name.toLowerCase()) {
            case "nn_table":
                nn_table = value;
                break;
            case "tmm_table":
                tmm_table = value;
                break;
            case "imm_table":
                imm_table = value;
                break;
            case "de_table":
                de_table = value;
                break;
            case "dnac1":
                dnac1 = Double.parseDouble(value);
                break;
            case "dnac2":
                dnac2 = Double.parseDouble(value);
            case "na":
                na = Double.parseDouble(value);
                break;
            case "k":
                k = Double.parseDouble(value);
                break;
            case "tris":
                tris = Double.parseDouble(value);
                break;
            case "mg":
                mg = Double.parseDouble(value);
                break;
            case "dntps":
                dntps = Double.parseDouble(value);
                break;
            case "saltcorr":
                saltcorr = Double.parseDouble(value);
                break;
        }
        //the script wants the primer to have a higher concentration than the template
        if (dnac2 > dnac1) {
            Double c = dnac2;
            dnac2 = dnac1;
            dnac1 = c;
        }
    }

    public static String getCmdArgs() {
        //the script wants the primer to have a higher concentration than the template
        if (dnac2 > dnac1) {
            Double c = dnac2;
            dnac2 = dnac1;
            dnac1 = c;
        }
        String s = "";
        String[] spar = {nn_table, tmm_table, imm_table, de_table};
        String[] sparnames = {"nn_table", "tmm_table", "imm_table", "de_table"};
        Double[] dpar = {dnac1, dnac2, na, k, tris, mg, dntps, saltcorr};
        String[] dparnames = {"dnac1", "dnac2", "Na", "K", "Tris", "Mg", "dNTPs", "saltcorr"};
        for (int i = 0; i < spar.length; i++) {
            if (spar[i] != null) {
                s = s + " " + sparnames[i] + "=" + spar[i];
            }
        }
        for (int i = 0; i < dpar.length; i++) {
            if (dpar[i] != null) {
                s = s + " " + dparnames[i] + "=" + dpar[i];
            }
        }
        return s;
    }
}
