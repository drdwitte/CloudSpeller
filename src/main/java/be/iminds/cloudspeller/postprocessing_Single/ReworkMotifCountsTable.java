package be.iminds.cloudspeller.postprocessing_Single;


import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by ddewitte on 03.01.15.
 */
public class ReworkMotifCountsTable {


    private static final int [] t = {15,50,60,70,90,95};

    private static final int [] cThr = { 50,70,90};
    private static final int [] blsThr = {15,70,95};
    private static final int [] Fthr = {1,5,10,15,20,500};
    private static final int [] motifLength = {6,7,8,9,10,11,12};
    private static final int [] degeneracies = {1,2,4,8,16,32,64};


    //key format: K S T C F
    public static void main(String [] args) throws IOException{

        String prefix = "AFRandom";
        String filename = "/home/ddewitte/Desktop/"+prefix;
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String line;

        Map<String,Integer> content = new HashMap<>();
        while ((line = br.readLine()) != null){
            String [] spl = line.split("\t");

            String key = spl[0];
            Integer value = Integer.parseInt(spl[1]);

            content.put(key,value);
        }

        br.close();


        //output: C, F, T k6S1 k6S4

        String outputFilename = "/home/ddewitte/Desktop/"+prefix+"Reworked.txt";
        BufferedWriter wr = new BufferedWriter(new FileWriter(outputFilename));


        for (int C : cThr){
            for (int F: Fthr){
                for (int T: blsThr){



                    wr.append(prefix+"_"+C+"_"+F+"_"+T+",");
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

        int c=0;
        while (c >= 0){
            c++;


        }

        for (int i=0; i<10; i++){
            System.out.println(c++);
        }
    }
}
