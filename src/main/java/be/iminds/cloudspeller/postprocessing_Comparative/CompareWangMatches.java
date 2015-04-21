package be.iminds.cloudspeller.postprocessing_Comparative;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import be.iminds.cloudspeller.toolbox.LineIterator;



public class CompareWangMatches {

	
	public static Map<String,Double> readDistanceMap(String filename) throws IOException{
		
		Map<String,Double> distanceMap = new HashMap<String,Double>();
		
		LineIterator it = new LineIterator(filename);
		
		while (it.hasNext()){
			String line = it.next();
			Scanner scanner = GeneralToolbox.generateScanner(line);
			
			String motif = scanner.next();
			String matchDistance = scanner.next();
			double distance = Double.parseDouble(matchDistance.split("_")[1]);
			
			if (distanceMap.get(motif)==null){
				distanceMap.put(motif,distance);
			}
			
		}
		
		
		return distanceMap;
	}
	
	public static void main (String [] args) throws IOException{
			
		String filename1 = args[0];
		String filename2 = args[1];
		
		String output = null;
		
		Map<String,Double> distanceMap1 = readDistanceMap(filename1);
		Map<String,Double> distanceMap2 = readDistanceMap(filename2);

		
		
		String filename3 = args[2];
		
		Set<String> wangMotifs = readWangMotifs(filename3);
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(output));
		out.write("Distance comparison: D_"+filename2+" - D_"+filename2+"\n");
		for (String wang : wangMotifs){
			
			double d1 = distanceMap1.get(wang);
			double d2 = distanceMap2.get(wang);
			
			out.write(wang+"\t"+(d2-d1)+"\n");
		}
		
		out.close();
		
	}

	private static Set<String> readWangMotifs(String filename) throws IOException {

		Set<String> motifs = new HashSet<String>();
		LineIterator it = new LineIterator(filename);

		while (it.hasNext()){
			String line = it.next();
			motifs.add(line.trim());
		}
		return motifs;
		
		
	}
}
