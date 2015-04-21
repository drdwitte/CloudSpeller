package be.iminds.cloudspeller.output;

public class Probability {

	private double prob;
	private static final int toPercent=100;
	
	public Probability(int occ, double backgroundOccs){
		if (occ>backgroundOccs){
			prob=toPercent*(occ - backgroundOccs)/occ;
		} else {
			prob=0.0;
		}
	}
	
	public Probability(double prob) {
		this.prob=prob;
	}

	public String toString(){
		return ""+prob;
	}
	
	public double getValue(){
		return prob;
	}
	
	
}
