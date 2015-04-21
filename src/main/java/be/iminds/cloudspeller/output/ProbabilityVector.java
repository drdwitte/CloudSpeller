package be.iminds.cloudspeller.output;

import be.iminds.cloudspeller.motifmodels.FreqVec;

public class ProbabilityVector {

	Probability [] probs;
	
	public ProbabilityVector(FreqVec freqVec, double[] backgroundOccs){
		this();
		for (int i=0; i<probs.length; i++){
			probs[i]=new Probability(freqVec.getFreq(i), backgroundOccs[i]);
		}
	}

	public ProbabilityVector() {
		probs = new Probability[FreqVec.getNumberOfIntervals()];
	}

	public Probability get(int i) {
		return probs[i];
	}
	
	public void setProb(int i, double prob){
		probs[i]=new Probability(prob);
	}

}
