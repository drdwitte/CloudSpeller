package be.iminds.cloudspeller.postprocessing_Single;

public class PoisssonEstimator {

	
	/**Calculates an approximation of the Poisson probability.
	 * @param mean - lambda, the average number of occurrences
	 * @param observed - the actual number of occurences observed
	 * @return ln(Poisson probability) - the natural log of the Poisson probability.
	 */
	public static double poissonProbabilityApproximation (double mean, int observed) {
		double k = observed;
		double a = k * Math.log(mean);
		double b = -mean;
		return a + b - factorialApproximation(k);
	}
	 
	/**Srinivasa Ramanujan ln(n!) factorial estimation.
	 * Good for larger values of n.
	 * @return ln(n!)
	 */
	public static double factorialApproximation(double n) {
		if (n < 2.) return 0;
		double a = n * Math.log(n) - n;
		double b = Math.log(n * (1. + 4. * n * (1. + 2. * n))) / 6.;
		return a + b + Math.log(Math.PI) / 2.;
	}
	
	/*public static double poissonCumulative(double mean, int observed){
		
		double p = 0.0; 
		double sum =0.0;
		int o=observed;
		double newP = poissonProbabilityApproximation(mean, o++);
		
		do {
			p=newP;
			sum+=p;
			newP = poissonProbabilityApproximation(mean, o++);
		}  while (newP > sum/100);
		  
		  
		return sum;  
		
		
		
	}*/
	
	public static double poissonApproxWiki(double mean, double observed){
		
		double f1 = Math.exp(-mean);
		double f2 = Math.pow(Math.E*mean,observed);
		double f3 = Math.pow(observed,observed);
		
		return f1*f2/f3;
		
		
	}
	
	
	public static void main (String [] args){
		
		int [] F = {1,2,3,4,5,10,20,50,100,250,500,1000,2000,5000};
		

		System.out.println("90% confidence threshold");

		
		for (int f : F){
			int observed = 10*f;
			double pValue = poissonApproxWiki(f,observed);
			System.out.println(f+"\t"+observed+"\t"+pValue);
			
		}
		
		
		System.out.println("50% confidence threshold");

		for (int f : F){
			int observed = 2*f;
			double pValue = poissonApproxWiki(f,observed);
			System.out.println(f+"\t"+observed+"\t"+pValue);
			
		}
		
	}
}
