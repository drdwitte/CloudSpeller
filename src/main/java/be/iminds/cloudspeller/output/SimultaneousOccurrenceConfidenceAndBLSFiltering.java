package be.iminds.cloudspeller.output;


import be.iminds.cloudspeller.phylogenetics.BLS;

public class SimultaneousOccurrenceConfidenceAndBLSFiltering extends SimultaneousOccurrenceAndConfidenceFiltering {


	public SimultaneousOccurrenceConfidenceAndBLSFiltering(
			int confidenceCutoff, int familyOccurrenceCutoff, int blsCutoff) {
		super(confidenceCutoff,familyOccurrenceCutoff);
		setBLSCutoff(blsCutoff);
		
	}
	
	//GETTERS
	public int getBLSCutoff(){
		return blsCutoff;
	}
	
	//SETTERS
	public void setBLSCutoff(int blsCutoff){
		this.blsCutoff=blsCutoff;
		calculateMinBLSID();
	}

	private void calculateMinBLSID() {
		
		int [] t = BLS.getBLSThresholds();
		
		if (blsCutoff>t[t.length-1]){ //maxID is length -1 -> BLS95 = BLS100
			minBLSId=t.length-1;
		}
		
		int i=0;
		
		while (blsCutoff > t[i]){ //BLS 20 maps to BLS 50 (in case of 15,50,60,70,90,95)
			i++;
		}
	
		minBLSId=i;
	}
	
	
}
