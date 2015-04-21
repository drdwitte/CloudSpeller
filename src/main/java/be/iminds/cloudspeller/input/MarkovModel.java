package be.iminds.cloudspeller.input;

import org.apache.commons.lang.NotImplementedException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by ddewitte on 21.04.15.
 */
public class MarkovModel {

    private static final String chars = "ACGT";
    private static final double smallProb = 1e-06;

    private Map<String, Integer> freqKmers;
    private Map<String, Integer> freqKm1Mers;


    private int numKmers=0;
    private int numKm1Mers=0;

    private int order;

    public MarkovModel(int order){

        if (order<1)
        throw new NotImplementedException("Order > 1");

        this.order=order;
        this.freqKmers = generateFreqMap(order+1);
        this.freqKm1Mers = generateFreqMap(order);

    }

    private static Map<String, Integer> generateFreqMap(int length) {

        List<String>  kmers = generateAllPossibleKMers(length);
        Map<String, Integer> kmerMap = new HashMap<>();

        for (String s : kmers){
            kmerMap.put(s,0);
        }

        return kmerMap;

    }

    private static List<String> generateAllPossibleKMers(int length){
        List<String> temp = new ArrayList<>();
        List<String> result = new ArrayList<>();

        temp.add("");

        for (int i=0; i<length; i++){
            for (String s : temp){
                for (Character c : chars.toCharArray()){
                    result.add(s + c);
                }
            }

            temp.clear();
            temp.addAll(result);
            result.clear();
        }

        return temp;


    }

    private void addKMer(String kmer){
        if (!nonAlph(kmer)){
            freqKmers.put(kmer,freqKmers.get(kmer)+1);
            numKmers++;
        }
    }
    private void addKm1Mer(String km1mer){
        if (!nonAlph(km1mer)){
            freqKm1Mers.put(km1mer,freqKm1Mers.get(km1mer)+1);
            numKm1Mers++;
        }
    }


    public void feed(String seq){

        int numberOfKmers = seq.length()-order;

        for (int i=0; i<numberOfKmers; i++){

            String kMer = seq.substring(i, i + order + 1);
            String kM1Mmer = seq.substring(i,i+order);

            addKm1Mer(kM1Mmer);
            addKMer(kMer);
        }

        //String kM1Mmer = seq.substring(numberOfKmers, numberOfKmers + order);
        //addKm1Mer(kM1Mmer);

    }

    public static boolean nonAlph(String s){

        for (Character c : s.toCharArray()){
            if (chars.indexOf(c)<0)
                return true;
        }
        return false;

    }

    public static void pr(String s){
        System.out.println(s);
    }


    public String generateSequence(int length){

        StringBuilder sb = new StringBuilder("");

        int prefixLength = Math.min(length,order);

        sb.append(generateRandomString(prefixLength));

        for (int i=order; i<length; i++){
            sb.append(generateCharacterForPrefix(sb.substring(i-order,i)));

        }
        return sb.toString();

    }

    public char generateCharacterForPrefix(String prefix){

        double [] probs = generateCharProbs(prefix);
        double p = Math.random();

        double sum = 0;
        for (int i=0; i<chars.length()-1; i++){
            sum+=probs[i];

            if (p<=sum){
                return chars.charAt(i);
            }
        }


        return chars.charAt(chars.length()-1);

    }

    private double[] generateCharProbs(String prefix) {

        double [] probs = new double[chars.length()];

        for (int i=0; i<probs.length; i++){



            double pABC = (1.0*freqKmers.get(prefix+chars.charAt(i)))
                    / numKmers;

            double pABC_renorm = (pABC + smallProb) / (1.0 + freqKmers.size()*smallProb);

            double pAB = (1.0*freqKm1Mers.get(prefix))
                    / numKm1Mers;

            double pAB_renorm = (pAB + smallProb) / (1.0 + freqKm1Mers.size()*smallProb);

            //pr(pABC + " -> " + pABC_renorm+"");
            //pr(pAB + " -> " + pAB_renorm+"");

            probs[i]= pABC_renorm / pAB_renorm;

            //pr(prefix+" -> "+chars.charAt(i)+" : " +probs[i]);

        }

        return probs;


    }

    public static char generateRandomCharacter(){
        int randomIndex = (int)(Math.random()*chars.length());
        return chars.charAt(randomIndex);
    }

    public static String generateRandomString(int length){
        StringBuilder sb = new StringBuilder("");
        for (int i=0; i<length; i++){
            sb.append(generateRandomCharacter());
        }
        return sb.toString();
    }


    public String toString(){

        return freqKm1Mers.toString() +"\n" + numKm1Mers
                + "\n" + freqKmers.toString() +"\n" + numKmers;

    }


    public static void main (String [] args){

        MarkovModel model = new MarkovModel(1);

        model.feed("ACGTACGTACGTGGAAAAAAAAAAAAAAAAAAAAAAAAAAAATTCTGAGAGAGAGAGACGY");

        System.out.println(model);

        int numSeq = 10000;
        long start = System.currentTimeMillis();
        List<String> l = new ArrayList<>();
        for (int i=0; i<numSeq; i++) {
            l.add(model.generateSequence(2000));
        }
        pr("done in "+(System.currentTimeMillis()-start));


        //5s voor 10000 sequenties =>


    }

}
