package be.iminds.cloudspeller.postprocessing_Single;


import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ddewitte on 03.01.15.
 */
public class ReworkMotifCountsTable {

    //use of logger: http://www.vogella.com/tutorials/Logging/article.html
    private final static Logger LOGGER = Logger.getLogger(ReworkMotifCountsTable.class.getName());

    private static final int [] t = {15,50,60,70,90,95};

    private static final int [] cThr = { 50,60,70,80,90,95};
    private static final int [] blsThr = {15,70,95};
    private static final int [] Fthr = {1,5,10,15,20,50,500};
    private static final int [] motifLength = {6,7,8,9,10,11,12};
    private static final int [] degeneracies = {1,2,4,8,16,32,64};


    //key format: K S T C F
    public static void main(String [] args) throws IOException{

        LOGGER.setLevel(Level.INFO);

        int alg = 1;
        int markov = 0;

        String path = "/home/ddewitte/FULL_PROFESSIONAL_BACKUP/Research/Bioinformatics_ALLSOURCES/" +
                "Bioinformatics@LaptopBureaublad/FDRTablesCFTKS_latestJuly2015/";

        String [] algs = {"AF/freqAF", "AB/freqAB" };
        String [] data = { "Real",  "Markov1",  "Markov2",  "Markov3" };
        String filenameLocal = "/part-r-00000";
        String filename = path + algs[alg] + data[markov] + filenameLocal;


        LOGGER.info("Reading: "+filename);

        BufferedReader br = new BufferedReader(new FileReader(filename));

        String line;

        Map<String,Long> content = new HashMap<>();
        while ((line = br.readLine()) != null){
            String [] spl = line.split("\t");

            String key = spl[0];
            Long value = Long.parseLong(spl[1]);

            if (value < 0){
                LOGGER.warning("Negative frequency corrected by add 2*Integer.MAX+1");
                LOGGER.warning("FROM: "+value);
                value=(long)Integer.MAX_VALUE*2+1+value;
                LOGGER.warning("TO: "+value);

            }

            content.put(key,value);
        }

        br.close();
        LOGGER.info("Finished Reading: "+filename);


        //output: C, F, T k6S1 k6S4

        String type = (algs[alg] + data[markov]).substring(7);
        String outputFilename = path + "/" + type + "Reworked.csv";


        LOGGER.info("Creating: "+outputFilename);

        BufferedWriter wr = new BufferedWriter(new FileWriter(outputFilename));

        String csvheader = "name,C,F,T," +
                "k6s1,k6s2,k6s4,k6s8,k6s16,k6s32,k6s64," +
                "k7s1,k7s2,k7s4,k7s8,k7s16,k7s32,k7s64," +
                "k8s1,k8s2,k8s4,k8s8,k8s16,k8s32,k8s64," +
                "k9s1,k9s2,k9s4,k9s8,k9s16,k9s32,k9s64," +
                "k10s1,k10s2,k10s4,k10s8,k10s16,k10s32,k10s64," +
                "k11s1,k11s2,k11s4,k11s8,k11s16,k11s32,k11s64," +
                "k12s1,k12s2,k12s4,k12s8,k12s16,k12s32,k12s64";

        wr.append(csvheader);
        wr.newLine();

        for (int C : cThr){
            for (int F: Fthr){
                for (int T: blsThr){



                    wr.append(type+"_"+C+"_"+F+"_"+T+",");
                    wr.append(C+",");
                    wr.append(F+",");
                    wr.append(T+",");

                    for (int k: motifLength){
                        for (int s: degeneracies){

                            String key = k+"_"+s+"_"+T+"_"+C+"_"+F;
                            wr.append(content.get(key)+",");
                        }
                    }

                    wr.newLine();

                }
            }
        }

        wr.close();

        LOGGER.info("Finished writing: "+outputFilename);






    }
}
