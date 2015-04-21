package be.iminds.cloudspeller.motifalgorithms;

import java.util.ArrayList;

import org.apache.commons.lang.NotImplementedException;

import be.iminds.cloudspeller.driver.MotifExtractor;

import be.iminds.cloudspeller.motifmodels.MotifFactory;

import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculator;

import be.iminds.cloudspeller.indexing.IndexStructureFactory;

import be.iminds.cloudspeller.input.Sequence;

public class DeNovoInexactDiscoveryAlgorithm implements DiscoveryAlgorithm {

	public DeNovoInexactDiscoveryAlgorithm(ArrayList<Sequence> arrayList) {
		throw new NotImplementedException();
	}

	@Override
	public void setDataStructure(IndexStructureFactory indexStructureFactory) {
		throw new NotImplementedException();
	}

	@Override
	public void setSearchSpace(MotifSearchSpace motifSearchSpace) {
		throw new NotImplementedException();
	}


	@Override
	public void setConservationScoreCalculator(ConservationScoreCalculator c) {
		throw new NotImplementedException();
	}


	@Override
	public void runDiscovery(MotifFactory factory) {
		throw new NotImplementedException();
		
	}

	@Override
	public void setMotifExtractor(MotifExtractor extractor) {
		throw new NotImplementedException();
	}






}
