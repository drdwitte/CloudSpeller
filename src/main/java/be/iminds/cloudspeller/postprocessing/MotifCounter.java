package be.iminds.cloudspeller.postprocessing;

import java.io.IOException;
import java.util.ArrayList;

import be.iminds.cloudspeller.output.MotifBLSRestrictions;

import be.iminds.cloudspeller.phylogenetics.BLS;

public class MotifCounter {

	private ArrayList<String> dataset;
	
	public MotifCounter(ArrayList<String> dataset){
		int [] thr = {15,50,60,70,90,95};
		BLS.initializeBLSConstants(thr);
		this.dataset=dataset;
	}
	
	public int [] count(MotifBLSRestrictions restrictions){
		
		int [] counts = new int [BLS.getBLSThresholds().length];  
		
		
		
		for (String file : dataset){
		
			try {
				
				BLSConfidenceGraphIterator iterator = new BLSConfidenceGraphIterator(file);
				
				while (iterator.hasNext()){
					boolean [] success = restrictions.checkRestrictions(iterator.next());
					
					for (int i=0; i<counts.length; i++){
						if (success[i]){
							counts[i]++;
						}
					}
				}
			
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return counts;
		
		
		
	}
}
