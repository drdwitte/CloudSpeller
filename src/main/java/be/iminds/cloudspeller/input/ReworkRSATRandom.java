package be.iminds.cloudspeller.input;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ddewitte on 05.06.15.
 */
public class ReworkRSATRandom {


    /**
     * Map the shuffled sequences back to a tfa format which can then be further pushed into gene families using
     * input.DatasetPreparator
     * @param filenameOut
     * @param randSeqFilename
     * @param geneSeqLengthsFilename
     */
    public static void generateTFA(String filenameOut, String randSeqFilename, String geneSeqLengthsFilename) throws IOException {

        BufferedWriter wr = new BufferedWriter(new FileWriter(new File(filenameOut)));
        List<String> geneNames = extractGeneNames(geneSeqLengthsFilename);
        List<String> sequences = extractSequences(randSeqFilename);

        System.out.println(geneNames.size() + " = " + sequences.size() + "!?");

        for (int i=0; i<geneNames.size(); i++){

            wr.write(">"+geneNames.get(i));
            wr.newLine();
            wr.write(sequences.get(i));
            wr.newLine();
        }


        wr.close();

    }

    private static List<String> extractSequences(String randSeqFilename) throws IOException  {
        List<String> l = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(new File(randSeqFilename)));

        StringBuilder tempSeq = new StringBuilder();

        String line;
        while ( (line = reader.readLine()) !=null){

            if (line.startsWith(">")){ //finish builder
                if (tempSeq.length()>0){
                    l.add(tempSeq.toString().toUpperCase());
                    tempSeq.setLength(0);
                }
            } else {
                tempSeq.append(line);
            }

        }

        if (tempSeq.length()>0){
            l.add(tempSeq.toString());
            tempSeq.setLength(0);
        }


        reader.close();
        return l;
    }

    private static List<String> extractGeneNames(String geneSeqLengthsFilename) throws IOException {
        List<String> l = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(new File(geneSeqLengthsFilename)));


        String line;

        boolean startReading = false;
        while ( (line = reader.readLine()) !=null){

            if (line.startsWith("#") || line.startsWith(";")){
                //ignore
            } else {
                String gene = line.split("\t")[0];
                l.add(gene);
            }

        }



        reader.close();
        return l;
    }

    public static void main (String [] args) throws IOException {


        String dir = "/home/ddewitte/Desktop/bioinformaticsPHD/OrigFormaatDataset/RSAT_shuffled";

        String  markovDir = "/Markov";



        int numMarkovOrders = 7;
        String [] species = {"bdi" , "zma", "sbi", "osa"};



        for (String sp : species){
            for (int m=0; m<numMarkovOrders; m++){


                String fullDirPath =  dir + markovDir + m + "/";

                String geneSeqLengthsFilename = fullDirPath + sp +".lengths";
                String randSeqFilename = fullDirPath + sp + "_markov" + m +".fasta";
                ;
                String filenameOut = fullDirPath + sp + "_markov" + m + ".tfa";
                ;
                ;

                System.out.println("Preparing: "+ filenameOut);
                generateTFA(filenameOut, randSeqFilename, geneSeqLengthsFilename);

            }
        }




        System.out.println("Done");
    }

}
