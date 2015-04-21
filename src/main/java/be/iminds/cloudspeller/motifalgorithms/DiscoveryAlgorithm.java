package be.iminds.cloudspeller.motifalgorithms;

import be.iminds.cloudspeller.driver.MotifExtractor;
import be.iminds.cloudspeller.motifmodels.MotifFactory;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculator;
import be.iminds.cloudspeller.indexing.IndexStructureFactory;

public interface DiscoveryAlgorithm extends MotifAlgorithm {
	
		public void runDiscovery(MotifFactory factory);
		
		public void setSearchSpace(MotifSearchSpace motifSearchSpace);
		public void setDataStructure(IndexStructureFactory indexStructureFactory);
		public void setConservationScoreCalculator(
				ConservationScoreCalculator createCalculator);
		public void setMotifExtractor(MotifExtractor extractor);

}

