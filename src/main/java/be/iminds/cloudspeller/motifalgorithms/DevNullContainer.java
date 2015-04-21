package be.iminds.cloudspeller.motifalgorithms;

import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.driver.MotifExtractor;

public class DevNullContainer implements MotifExtractor { 

	int numberOfMotifs;

	public DevNullContainer(){
		this.numberOfMotifs=0;
	}
	
	@Override
	public void add(Motif motif, ConservationScore score) {
		numberOfMotifs++;
	}
	
	@Override
	public int getNumberOfMotifsExtracted() {
		return numberOfMotifs;
	}
	
	@Override
	public void reset() {
		numberOfMotifs=0;
	}

	@Override
	public void close() {
		reset();
	}
	
	
}
