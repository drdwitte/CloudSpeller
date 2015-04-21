package be.iminds.cloudspeller.phylogenetics;

import org.apache.commons.lang.NotImplementedException;

import be.iminds.cloudspeller.motifmodels.FreqVec;

public class SC implements ConservationScore {

	//private static int numberOfOrganisms;
	//private static SortedMap<String> orgs;

	
	
	//GETTERS
	/*public static void initializeSCConstants(int minNumberOfSpecies, int nOrg) {
		numberOfOrganisms=nOrg;
		
	}*/
	
	//SETTERS
	
	//METHODS
	
	@Override
	public FreqVec createFrequencyVector() {
		throw new NotImplementedException();
	}

	@Override
	public int compareTo(ConservationScore o) {
		throw new NotImplementedException();
	}
	
	/**
	 * Recursive method to calculate combinatorial C^k_n = n! / k! (n-k)!
	 * C^k_n = C^{k-1}_{n-1}
	 * @param k
	 * @param n
	 * @return
	 */
	/*private static int combinatorial(int k, int n){
		int denominator=1;
		int numerator=1;
		for (int i=0; i<k-1; i++){
			denominator*=(n-i);
			numerator*=(k-i);
		}
		return denominator/numerator;
	}*/

	


}
