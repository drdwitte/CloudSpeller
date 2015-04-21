package be.iminds.cloudspeller.motifalgorithms;

import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.driver.ABMotifExtractor;

public class ABMotifContainer implements ABMotifExtractor {

	int numberOfMotifsExtracted=0;
	private Map<Motif,ConservationScore> motifMap;


	public ABMotifContainer(){
		motifMap = new HashMap<Motif,ConservationScore>();
	}

	@Override
	/**
	 * NOTE: pay attention to reverse motif must also be added!!
	 */
	public void add(Motif motif, ConservationScore score) {
			
		ConservationScore oldScore = motifMap.get(motif);
				
		if (oldScore == null || oldScore.compareTo(score)<0){
			motifMap.put(motif,score);

		} else {}
		
	}

	@Override
	public void close() {
		reset();
	}

	@Override
	public int getNumberOfMotifsExtracted() {
		return numberOfMotifsExtracted;
	}

	@Override
	/**
	 * Clear motif map en emit motifs (after maptask has completed
	 */
	public void reset() {
		motifMap.clear();
	}

	public Map<Motif, ConservationScore> getMotifMap() {
		return motifMap;
	}


}



