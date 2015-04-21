package be.iminds.cloudspeller.phylogenetics;

import be.iminds.cloudspeller.motifmodels.FreqVec;

public class BLS implements ConservationScore {
	
	public static int MIN=10;
	private static int[] blsThresholds = {10,30,50,70,90}; //default thresholds
	
	private int bls;
	
	@Override
	public FreqVec createFrequencyVector() {
		int one=1;
		FreqVec v= new FreqVec();
		for (int i=0; i<blsThresholds.length; i++){
			if (bls>=blsThresholds[i]){
				v.setFreq(i,one);
			} else {
				return v;
			}
		}
		return v;
	}
	
	public static void initializeBLSConstants(int min, int width, int numberOfIntervals){
		blsThresholds = new int [numberOfIntervals];
		
		for (int i=0; i<blsThresholds.length; i++){
			blsThresholds[i]=min+i*width;
		}
		MIN = blsThresholds[0];
		
	}
	
	/**
	 * Take a deep copy of the bls thresholds in the argument
	 * @param thresholds
	 */
	public static void initializeBLSConstants(int [] thresholds){
		blsThresholds = new int [thresholds.length];
		for (int i=0; i<blsThresholds.length; i++){
			blsThresholds[i]=thresholds[i];
		}
		MIN = blsThresholds[0];
	}
	
	
	public static int getNumberOfIntervals(){
		return blsThresholds.length;
	}
	
	
	public BLS(int score){
		this.bls=score;
	}

	public String toString(){
		return ""+bls;
	}

	@Override
	public int compareTo(ConservationScore score) {

		if (!(score instanceof BLS)){
			throw new ClassCastException("Invalid conservation score comparison");
		}
		BLS BLSScore = (BLS) score;
		return (this.bls-BLSScore.bls);
	}
			
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof BLS)){
			throw new ClassCastException("Invalid conservation score comparison");
		}
		
		BLS BLSScore = (BLS) obj;
		return (this.bls==BLSScore.bls);
	}

	public static int[] getBLSThresholds() {
		return blsThresholds;
	}

	public double doubleValue() {
		return 1.0*bls;
	}
	
	public int intValue(){
		return bls;
	}

}



	

