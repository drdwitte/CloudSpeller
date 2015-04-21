package be.iminds.cloudspeller.postprocessing;

import be.iminds.cloudspeller.indexing.PatternBLSPair;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.MotifBLSRestrictions;
import be.iminds.cloudspeller.phylogenetics.BLS;

public class OutputComparator {

	private String dataset1;
	private String dataset2;

	public OutputComparator(String dataset1, String dataset2){
		this.dataset1=dataset1;
		this.dataset2=dataset2;
	}
	
	public Set<String> calculateMotifIntersection(ConfidenceGraphRestrictions restrictions){
		Set<String> intersection = new HashSet<String>();
		
		Set<String> motifsD1 = new HashSet<String>();
		
		try {
			BLSConfidenceGraphIterator iterator = new BLSConfidenceGraphIterator(dataset1);
			
			while (iterator.hasNext()){
				BLSConfidenceGraph graph = iterator.next();
				if (restrictions.checkRestrictions(graph)){
					motifsD1.add(graph.getMotif().toString());
				}
			}
		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			BLSConfidenceGraphIterator iterator = new BLSConfidenceGraphIterator(dataset2);
			
			while (iterator.hasNext()){
				String motifD2 = iterator.next().getMotif().toString();
				
				if (motifsD1.contains(motifD2)){
					intersection.add(motifD2);
				}
			}
		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return intersection;
	}
	
	public Set<PatternBLSPair> calculateMotifBLSIntersection(MotifBLSRestrictions restrictions){
		Set<PatternBLSPair> intersection = new HashSet<PatternBLSPair>();
		Set<PatternBLSPair> motifBLSD1 = new HashSet<PatternBLSPair>();
		
		try {
			BLSConfidenceGraphIterator iterator = new BLSConfidenceGraphIterator(dataset1);
			
			while (iterator.hasNext()){
				BLSConfidenceGraph graph = iterator.next(); 
				boolean [] success = restrictions.checkRestrictions(graph);
				
				for (int i=0; i<success.length; i++){
					if (success[i]){
						motifBLSD1.add(new PatternBLSPair(graph.getMotif().toString(),BLS.getBLSThresholds()[i]));	
					}
				}
			}
		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			BLSConfidenceGraphIterator iterator = new BLSConfidenceGraphIterator(dataset2);
			
			while (iterator.hasNext()){
				
				BLSConfidenceGraph graph = iterator.next(); 
				boolean [] success = restrictions.checkRestrictions(graph);
				
				for (int i=0; i<success.length; i++){
					if (success[i]){
						PatternBLSPair pbls = new PatternBLSPair(graph.getMotif().toString(), BLS.getBLSThresholds()[i]);
						if (motifBLSD1.contains(pbls)){
							intersection.add(pbls);
						}	
					}
				}	
			}
		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		
		
		
		
		
		
		
		
		return intersection;
	}
	
}
