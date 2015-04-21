package be.iminds.cloudspeller.phylogenetics;



import be.iminds.cloudspeller.motifmodels.FreqVec;

public interface ConservationScore extends Comparable<ConservationScore> {
	
	FreqVec createFrequencyVector();
}
