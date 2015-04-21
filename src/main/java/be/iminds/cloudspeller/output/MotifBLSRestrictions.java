package be.iminds.cloudspeller.output;


import be.iminds.cloudspeller.phylogenetics.BLS;

import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import be.iminds.cloudspeller.motifmodels.FreqVec;

public class MotifBLSRestrictions {

	protected int famCutoff;
	protected double probCutoff;
	
	protected int blsCutoff;
	protected int minBLSId;
	protected int [] F;
	protected double [] p;
		

	public MotifBLSRestrictions(int famCutoff, double probCutoff){
		this.famCutoff=famCutoff;
		this.probCutoff=probCutoff;
		this.blsCutoff=0;
		this.minBLSId=0;
		F = new int[BLS.getBLSThresholds().length];
		p = new double[BLS.getBLSThresholds().length];
	}
	
	public boolean [] checkRestrictions(BLSConfidenceGraph graph){
		boolean [] success = new boolean [FreqVec.getNumberOfIntervals()];		
		for (int i=minBLSId; i<success.length; i++){
			int fi = graph.getFreq(i);
			double pi = graph.getProbValue(i);
			
			success[i]= fi >= famCutoff && pi >=probCutoff;
		}
		
		return success;
		
	}

	/*public int getMaxFamilies(BLSConfidenceGraph graph) {
		int F=-1;
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			int fi = graph.getFreq(i);
			double pi = graph.getProbValue(i);
			
			if (fi >= famCutoff && pi >=probCutoff ){
				return fi;
			}
		}
		return F;
	}*/
	
	
	
	
	
	public int getMaxFamilies(String outputLine){
		
		GeneralToolbox.parseConfidenceGraphValues(outputLine, F, p);
		int Fmax=-1;
	
		for (int i=minBLSId; i<F.length; i++){

			if (F[i] >= famCutoff && p[i] >=probCutoff ){
				return F[i];
			}
		}
		return Fmax;
	}

}
