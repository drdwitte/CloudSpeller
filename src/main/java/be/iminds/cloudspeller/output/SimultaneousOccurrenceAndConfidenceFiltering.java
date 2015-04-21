package be.iminds.cloudspeller.output;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import be.iminds.cloudspeller.motifmodels.FreqVec;

public class SimultaneousOccurrenceAndConfidenceFiltering implements ConfidenceGraphRestrictions{

	protected double probCutoff;
	protected int famCutoff;
	
	protected int blsCutoff;
	protected int minBLSId;
	protected int [] F;
	protected double [] p;
	
	/**
	 * Constructor
	 * @param confidenceCutoff
	 * @param familyOccurrenceCutoff
	 */
	public SimultaneousOccurrenceAndConfidenceFiltering(double probCutoff,
			int famCutoff) {
		
		this.famCutoff=famCutoff;
		this.probCutoff=probCutoff;
		this.blsCutoff=0;
		this.minBLSId=0;
		F = new int[BLS.getBLSThresholds().length];
		p = new double[BLS.getBLSThresholds().length];
	}


	//GETTERS
	
	public double getConfidenceCutoff() {
		return probCutoff;
	}

	public int getFamilyOccurrenceCutoff() {
		return famCutoff;
	}
	
	//SETTERS
	
	public void setConfidenceCutoff(int probCutoff) {
		this.probCutoff = probCutoff;
	}

	public void setFamilyOccurrenceCutoff(int familyOccurrenceCutoff) {
		this.famCutoff = familyOccurrenceCutoff;
	}
	
	
	//METHODS
	
	@Override
	public boolean checkRestrictions(ConfidenceGraph graph) {
		for (int i=minBLSId; i<FreqVec.getNumberOfIntervals(); i++){
			if (graph.getProbValue(i)>=probCutoff){
				if (graph.getFreq(i)>=famCutoff)
					return true;
			}
		}
		return false;
	}


	@Override
	public boolean checkRestrictions(String strValue) {
		
		GeneralToolbox.parseConfidenceGraphValues(strValue, F, p);
		for (int i=minBLSId; i<FreqVec.getNumberOfIntervals(); i++){
			if (p[i]>=probCutoff){
				if (F[i]>=famCutoff)
					return true;
			}
		}
		return false;
		
	}

	@Override
	public boolean checkRestrictions(int [] F, double [] p){
		this.F = F;
		this.p = p;
		for (int i=minBLSId; i<FreqVec.getNumberOfIntervals(); i++){
			if (p[i]>=probCutoff){
				if (F[i]>=famCutoff)
					return true;
			}
		}
		return false;

	}










}
