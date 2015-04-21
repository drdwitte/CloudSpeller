package be.iminds.cloudspeller.indexing;


import be.iminds.cloudspeller.motifmodels.MotifFactory;

public interface IndexStructure extends PatternMatcher{
	
	public ISMonkey getExactISMonkey(MotifFactory motifFactory, int maxDegeneracy);
	public ISMonkey getInexactISMonkey();
}
