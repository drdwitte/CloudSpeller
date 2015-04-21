package be.iminds.cloudspeller.postprocessing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;


public class OutputReworker {

	private Alphabet alphabet = new IUPACAlphabet(IUPACType.TWOFOLDSANDN);
	
	
	private Map<String,ConfGraphCollector> generateOutputCollectors(String prefix){
		Map<String,ConfGraphCollector> outputCollectors = new HashMap<String,ConfGraphCollector>();
		
		for (Character c1 : alphabet){
			for (Character c2 : alphabet){
				for (Character c3 : alphabet){
					String fileID = ""+c1+c2+c3;
					outputCollectors.put(fileID,new ConfGraphCollector(1000,prefix+fileID+".txt"));
				}
			} 
		}
		
		return outputCollectors;
	}
	
	 
	
	//files will be labeled with motif prefix
	public void reworkOutputToContainMultiplePermGroups(String prefix, ArrayList<String> files) throws IOException{
		
		int fileIDLength=3;
		Map<String,ConfGraphCollector> outputCollectors =generateOutputCollectors(prefix);
		
		for (String file : files){
			BLSConfidenceGraphIterator iterator = new BLSConfidenceGraphIterator(file);
			while (iterator.hasNext()){
				BLSConfidenceGraph graph = iterator.next();
				String key = graph.getMotif().toString().substring(0,fileIDLength);
				outputCollectors.get(key).addGraph(graph);
			}
		}
		
		//flush final files
		for (Map.Entry<String,ConfGraphCollector> e : outputCollectors.entrySet()){
			e.getValue().flushCollector();
		}
	}
}




