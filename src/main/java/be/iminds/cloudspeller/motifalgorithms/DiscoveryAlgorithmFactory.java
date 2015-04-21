package be.iminds.cloudspeller.motifalgorithms;

import be.iminds.cloudspeller.input.Sequence;

import java.util.ArrayList;


public class DiscoveryAlgorithmFactory  {

	//benchmark paramters
	private static int maximalNumberOfMotifs=10000;
	private static BenchmarkType benchmarkType=BenchmarkType.EMITRANDOM;
	
	public static DiscoveryAlgorithm createAlgorithm(ArrayList<Sequence> seqs, MotifAlgType type) {
		switch(type){
			case EXACT:
				return new DeNovoExactDiscoveryAlgorithm(seqs);
			
			case EXACTAB:
				return new AlBasedExactDiscoveryAlgorithm(seqs);
				
			case INEXACT:
				return new DeNovoInexactDiscoveryAlgorithm(seqs);
			
			case FAKE:
				return new FakeDiscoveryAlgorithm(seqs,maximalNumberOfMotifs,benchmarkType);				
		
		}
		
		return null;
	}

	public static void setMaximalNumberOfMotifs(int maximalNumberOfMotifs) {
		DiscoveryAlgorithmFactory.maximalNumberOfMotifs = maximalNumberOfMotifs;
	}

	public static void setBenchmarkType(BenchmarkType benchmarkType) {
		DiscoveryAlgorithmFactory.benchmarkType = benchmarkType;
	}
	
}

