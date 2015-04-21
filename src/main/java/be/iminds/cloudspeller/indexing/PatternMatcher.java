package be.iminds.cloudspeller.indexing;

import java.util.List;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;

public interface PatternMatcher {

	public List<Suffix> matchExactPattern(IUPACMotif pattern);
}
