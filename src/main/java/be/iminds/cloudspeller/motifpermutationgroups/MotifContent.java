package be.iminds.cloudspeller.motifpermutationgroups;

import java.util.Set;

import be.iminds.cloudspeller.motifmodels.Motif;


public interface MotifContent {

	public Set<Motif> createPermutationGroup(int numberOfTries);

}
