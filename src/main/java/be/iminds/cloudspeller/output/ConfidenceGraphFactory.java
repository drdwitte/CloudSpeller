package be.iminds.cloudspeller.output;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;

public interface ConfidenceGraphFactory {

	ConfidenceGraph createGraph(Motif motif, FreqVec freqVec, ProbabilityVector probVec);
}
