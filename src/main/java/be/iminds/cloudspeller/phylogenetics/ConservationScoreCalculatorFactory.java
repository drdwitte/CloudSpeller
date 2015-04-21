package be.iminds.cloudspeller.phylogenetics;

import be.iminds.cloudspeller.input.GeneFamily;

public interface ConservationScoreCalculatorFactory {

	ConservationScoreCalculator createCalculator(GeneFamily gf);
}
