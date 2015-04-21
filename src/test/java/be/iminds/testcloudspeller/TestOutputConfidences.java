package be.iminds.testcloudspeller;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.*;
import org.junit.Test;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;

import be.iminds.cloudspeller.output.BLSConfGraphFactory;
import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.phylogenetics.BLS;

public class TestOutputConfidences {

	
	private static Map<Motif,BLSConfidenceGraph> output = new HashMap<Motif,BLSConfidenceGraph> ();
	private static Alphabet alphabet = new IUPACAlphabet(IUPACType.FULL);
	
	@Test
	public void testConfidences() throws IOException{
		
		FreqVec.setNumberOfIntervals(6);
		BLS.initializeBLSConstants(40,10,6);
		
		String [] filenames = {"ConfidencesAACKSTgroup.txt", "ConfidencesCCTTWYgroup.txt",
				"ConfidencesCGGGGNgroup.txt"};
		String [] permGroupIDs = {"AACKST", "CCTTWY", "CGGGGN"};
		
		
		
		for (int i=0; i<filenames.length; i++){
			testPermGroup(filenames[i],permGroupIDs[i]);
		}
	}
	
	private void testPermGroup(String filename, String permGroupID) throws IOException {
		
		BufferedReader in = new BufferedReader(new FileReader(filename));
		output.clear();
		new BLSConfGraphFactory().getGraphsFromBuffer(in, output);
		
		double backgr [] = generateBackground(getGroupSize(permGroupID));
		
		System.out.println(permGroupID);
		for (int i=0; i<backgr.length; i++){
			System.out.print(backgr[i]+" ");
		}
		System.out.println("");
		
		for (Map.Entry<Motif,BLSConfidenceGraph> e : output.entrySet()){
			BLSConfidenceGraph graph = e.getValue();
			
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				double fi =graph.getFreq(i); 
				double bi;
				double checkValue;
				
				
				if (fi>0.0){
					bi = backgr[i];
					checkValue =  100.0*(fi-bi)/fi;
					checkValue = (checkValue<0.0)?0.0:checkValue;
				} else {
					checkValue=0.0;
				}
				
				if (checkValue!=graph.getProbValue(i)){
					System.out.println(graph);
					System.out.println(checkValue);
					System.out.println(graph.getProbValue(i));
					System.out.println(i);
				}
				
				assertEquals(checkValue,graph.getProbValue(i),3);
			}
		}
		

	}
	
	private int getGroupSize(String s){
		int denominator = fac(s.length()); 
		double res = denominator;
		int [] charFreqs = getCharacterFrequencies(s);
		for (int i = 0; i<charFreqs .length; i++){
			res/=fac(charFreqs[i]);
		}
		return (int)res;
		
	}
	
	private int[] getCharacterFrequencies(String s) {
		int [] freqs = new int[alphabet.getAllChars().length()];
		String allChars = alphabet.getAllChars();
		
		for (int i=0; i<allChars.length(); i++){
			for (int j=0; j<s.length(); j++){
				if (s.charAt(j)==allChars.charAt(i)){
					freqs[i]++;
				}
			}
			
		}
		
		
		return freqs;
	}

	private int fac( int f ){
		int res = 1;
		for (int i=2; i<=f; i++){
			res*=i;
		}
		return res;
	}
	
	/**
	 * MEDIAN code from stackOverflow: 
	 * http://stackoverflow.com/questions/4191687/how-to-calculate-mean-median-mode-and-range-from-a-set-of-numbers
	 */
	// the array double[] m MUST BE SORTED
	public static double median(Double[] m) {
	    int middle = m.length/2;
	    if (m.length%2 == 1) {
	        return m[middle];
	    } else {
	        return (m[middle-1] + m[middle]) / 2.0;
	    }
	}
	
	/**
	 * NOTE: only works when maxGroupsize < 1000 which is the case for this test k=5-6 and degChars<=2
	 * => 6!/2! = 360 is the maximum number of control motifs
	 * @param maxGroupSize
	 * @return
	 */
	double [] generateBackground(int maxGroupSize){
		double [] backgr = new double [FreqVec.getNumberOfIntervals()]; 
		
		
		for (int i=0; i<backgr.length; i++){
			
			ArrayList<Double> freqList = new ArrayList<Double>();
			int c=0;
			for (Map.Entry<Motif,BLSConfidenceGraph> e : output.entrySet()){
				freqList.add((double)e.getValue().getFreq(i));
				c++;
			}
			
			//add zeroes for non occurring background motifs
			for (; c<maxGroupSize; c++){
				freqList.add(0.0);
			}
			
			
			Collections.sort(freqList);
			Double [] sortedArray = new Double[maxGroupSize]; 
			sortedArray = freqList.toArray(sortedArray); 
			backgr[i]=median(sortedArray);
		}
		
		
		return backgr;
	}
	
	
}
