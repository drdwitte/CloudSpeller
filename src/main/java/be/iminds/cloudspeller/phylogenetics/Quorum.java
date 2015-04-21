package be.iminds.cloudspeller.phylogenetics;

import org.apache.commons.lang.NotImplementedException;

import be.iminds.cloudspeller.motifmodels.FreqVec;

public class Quorum implements ConservationScore {

	@Override
	public FreqVec createFrequencyVector() {
		throw new NotImplementedException();
	}

	@Override
	public int compareTo(ConservationScore o) {
		throw new NotImplementedException();
	}
}
