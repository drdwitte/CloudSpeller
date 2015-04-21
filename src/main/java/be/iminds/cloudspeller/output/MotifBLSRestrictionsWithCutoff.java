package be.iminds.cloudspeller.output;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.phylogenetics.BLS;

public class MotifBLSRestrictionsWithCutoff extends MotifBLSRestrictions {

	public MotifBLSRestrictionsWithCutoff(int famCutoff, double probCutoff, int blsCutoff) {
		super(famCutoff, probCutoff);
		setBLSCutoff(blsCutoff);
		
		
	}
	
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
	
	@Override
	public boolean [] checkRestrictions(BLSConfidenceGraph graph){
		boolean [] success = new boolean [FreqVec.getNumberOfIntervals()];		
		for (int i=minBLSId; i<success.length; i++){
			int fi = graph.getFreq(i);
			double pi = graph.getProbValue(i);
			
			success[i]= fi >= famCutoff && pi >=probCutoff;
		}
		
		return success;
		
	}

	/*@Override
	public int getMaxFamilies(BLSConfidenceGraph graph) {
		int F=-1;
		for (int i=minBLSId; i<FreqVec.getNumberOfIntervals(); i++){
			int fi = graph.getFreq(i);
			double pi = graph.getProbValue(i);
			
			if (fi >= famCutoff && pi >=probCutoff ){
				return fi;
			}
		}
		return F;
	}*/
	


}
