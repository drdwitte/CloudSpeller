package be.iminds.cloudspeller.output;

import org.apache.avro.generic.GenericData;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by ddewitte on 10.07.15.
 */
public class KmerDistributionCorrelator {



    private static final String dirPrefix = "/home/ddewitte/Desktop/bioinformaticsPHD/FreqDistributionsRealandMarkov/";
    private static final String fileSuffix = "part-r-00000";

    //real data
    private static final String realStr = "_real";

    //markov 0 tem 4
    private static final String markovStr = "_markov";

    public static void main(String [] args) throws IOException {

        /*Map<String,Integer > test1 = new HashMap<>();
        Map<String,Integer > test2 = new HashMap<>();
        test1.put("a",2589);
        test1.put("c",2161);
        test2.put("a",2557);
        test2.put("c",2199);

        System.out.println(calculatePearsonCorrelation(test1,test2));
*/


        String [] orgs = {"bdi", "osa", "sbi", "zma"};

        int markov_min = 0;
        int markov_max = 4;

        int kmin = 1;
        int kmax = 7;




        for (String org: orgs){
            System.out.println("Pearson correlation results for "+org);
            System.out.println("k\tmarkov0\tmarkov1\tmarkov2\tmarkov3\tmarkov4");

            for (int k=kmin; k<=kmax; k++){

                StringBuilder sb = new StringBuilder();
                sb.append("k="+k);

                Map<String,Long> kmersReal = extractKmersFromFile(filenameReal(org),k);

                for (int m=markov_min; m<=markov_max; m++){

                    Map<String,Long> kmersMarkov = extractKmersFromFile(filenameMarkov(org,m),k);


                    double corr = calculatePearsonCorrelation(kmersReal,kmersMarkov);

                    sb.append("\t"+corr);

                }

                System.out.println(sb);


            }

            System.out.println("");

        }

        System.out.println("");
        System.out.println("");



        for (String org: orgs){
            System.out.println("Pearson correlation results for "+org + " Random permutation of real!");
            System.out.println("k\tmarkov0\tmarkov1\tmarkov2\tmarkov3\tmarkov4");

            for (int k=kmin; k<=kmax; k++){

                StringBuilder sb = new StringBuilder();
                sb.append("k="+k);

                Map<String,Long> kmersRealtemp = extractKmersFromFile(filenameReal(org),k);

                //pr("Voor shuffle: "+kmersRealtemp);

                Map<String,Long> kmersReal = shuffleKV(kmersRealtemp);

                //pr("Na shuffle: "+kmersReal);


                for (int m=markov_min; m<=markov_max; m++){

                    Map<String,Long> kmersMarkov = extractKmersFromFile(filenameMarkov(org,m),k);

                    //pr("Markov vector: "+kmersMarkov);

                    double corr = calculatePearsonCorrelation(kmersReal,kmersMarkov);

                    sb.append("\t"+corr);

                }

                System.out.println(sb);


            }

            System.out.println("");

        }



    }

    private static Map<String, Long> shuffleKV(Map<String, Long> kmers) {

        List<Long> values = new ArrayList<>();
        Map<String,Long> newKmers = new HashMap<>();

        for (String s : kmers.keySet()){

            values.add(kmers.get(s));

        }

        Collections.shuffle(values);

        int counter = 0;
        for (String s : kmers.keySet()){

            newKmers.put(s, values.get(counter++));
        }



        return newKmers;
    }

    private static String filenameReal(String org) {

        return dirPrefix + org + realStr + "/" + fileSuffix ;
    }

    private static String filenameMarkov(String org, int m) {

        return dirPrefix + org + markovStr + m +  "/" + fileSuffix ;

    }

    private static Map<String, Long> extractKmersFromFile(String filename, int k) throws IOException {

        Map<String, Long> kvmap = new HashMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line;
        while ( (line = reader.readLine()) != null){
            String [] parts = line.split("\t");
            String kmer = parts[0];

            if (kmer.length() == k){
                Long freq = Long.parseLong(parts[1]);
                kvmap.put(kmer,freq);
            }
        }

        return kvmap;

    }

    private static double Exy(Map<String, Long> v1, Map<String, Long> v2){

        double sum = 0.0;

        //pr("in"); pr(v1); pr(v2);


        int numDimensions = v1.size();

        if (v2.size()!=numDimensions){
            System.err.println("Some strings do not occur!");
        }

        for (Map.Entry<String,Long> e: v1.entrySet()){

            if (v2.get(e.getKey())==null){
                System.err.println("Komt niet voor!!!" + e.getKey());
            } else {

                sum += e.getValue() * v2.get(e.getKey());
                //pr(e.getValue()); pr(v2.get(e.getKey()));
                //pr(e.getValue() * v2.get(e.getKey()));
            }



        }

        return sum/numDimensions;
    }

    private static double Ex(Map<String, Long> v){

        double sum = 0.0;
        int numDimensions = v.size();

        for (Map.Entry<String,Long> e: v.entrySet()){

            sum+= e.getValue();

        }

        return sum/numDimensions;
    }



    private static double calculatePearsonCorrelation(Map<String, Long> v1, Map<String, Long> v2 ){

        double EXY = Exy(v1,v2);
        double EX = Ex(v1);
        double EY = Ex(v2);;
        double EX2 = Exy(v1,v1);
        double EY2 = Exy(v2,v2);

        /*pr("EXY: "+EXY);
        pr("EX: "+EX);
        pr("EY: "+EY);
        pr("EXX: "+EX2);
        pr("EYY: "+EY2);*/



        double nom = EXY - EX*EY;
        double denom1 = Math.sqrt(EX2 - EX*EX);
        double denom2 = Math.sqrt(EY2 - EY*EY);

        /*pr(nom);
        pr(EX2 - EX*EX);
        pr(EY2 - EY*EY);*/



        return nom / denom1 / denom2;

    }

    private static void pr(Object o){
        System.out.println(o.toString());
    }


}
