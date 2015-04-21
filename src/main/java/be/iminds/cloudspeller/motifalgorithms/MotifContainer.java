package be.iminds.cloudspeller.motifalgorithms;

import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.driver.MotifExtractor;
/**
 * Warning does not allow duplicate entries so no aggregation in add() method!!
 * @author ddewitte
 *
 */
public class MotifContainer implements MotifExtractor {

	Map<Motif,ConservationScore> motifMap = new HashMap<Motif,ConservationScore>();
	
	@Override
	public void add(Motif motif, ConservationScore score) {
		motifMap.put(motif,score);
	}

	@Override
	public int getNumberOfMotifsExtracted() {
		return motifMap.size();
	}

	public Map<Motif, ConservationScore> getMotifMap(){
		return motifMap;
	}

	@Override
	public void reset() {
		motifMap.clear();
	}

	@Override
	public void close() {
		reset();
	}
}
