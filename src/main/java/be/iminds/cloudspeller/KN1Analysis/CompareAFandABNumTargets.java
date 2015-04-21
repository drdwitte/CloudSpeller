package be.iminds.cloudspeller.KN1Analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import be.iminds.cloudspeller.toolbox.LineIterator;

public class CompareAFandABNumTargets {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		String AFfile = "/home/ddewitte/Bureaublad/CompareAFenABKN1/AFC90F0";
		String ABfile = "/home/ddewitte/Bureaublad/CompareAFenABKN1/ABC90F0";
		
		Map<String,Integer> motifNumTargetsAF = reworkFile(AFfile);
		Map<String,Integer> motifNumTargetsAB = reworkFile(ABfile);
	
		
		Set<String> allMotifs = new HashSet<String>();
		allMotifs.addAll(motifNumTargetsAB.keySet());
		allMotifs.addAll(motifNumTargetsAF.keySet());
		
		SortedMap<Integer,Integer> freqMap = new TreeMap<Integer,Integer>();
		Integer Nul = new Integer(0);
		for (String motif : allMotifs){
			
			Integer iAF=motifNumTargetsAF.get(motif);
			Integer iAB=motifNumTargetsAB.get(motif);
			
			if (iAF==null){
				iAF=Nul;
			}
			if (iAB==null){
				iAB=Nul;
			}
			
			int diff = iAF-iAB;
			
			Integer oldFreq = freqMap.get(diff);
			
			if (oldFreq==null){
				oldFreq=Nul;
			}
			
			freqMap.put(diff,oldFreq+1);
			
		}
		
		for (Map.Entry<Integer,Integer> e : freqMap.entrySet()){
			System.out.println(e.getKey()+"\t"+e.getValue());
		}
		
		
		
	}

	private static Map<String, Integer> reworkFile(String file) throws IOException {
		
		Map<String, Integer> m = new HashMap<String,Integer>();
		
		LineIterator iterator = new LineIterator(file);
		
		while (iterator.hasNext()){
			
			String line = iterator.next();
			
			if (line.length()==0)
				continue;
			
			String [] spl = line.split("\t");
			
			String motif = spl[0];
			int numMaizeTargets = Integer.parseInt(spl[6]);
			
			m.put(motif,numMaizeTargets);
			
		}
		
		return m;
	}

}
