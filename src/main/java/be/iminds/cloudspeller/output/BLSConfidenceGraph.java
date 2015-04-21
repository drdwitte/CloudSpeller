package be.iminds.cloudspeller.output;

import java.util.Scanner;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import org.apache.hadoop.io.Text;


import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;

public class BLSConfidenceGraph implements ConfidenceGraph {

	private static final char newline = '\n';
	private static final char tab = '\t';
	
	
	private Motif motif;
	private FreqVec freqVec;
	private ProbabilityVector probVec;
	
	public BLSConfidenceGraph(Motif motif, FreqVec freqVec, ProbabilityVector probVec){
		setMotif(motif);
		setFreqVec(freqVec);
		setProbVec(probVec);
	}
	
	public BLSConfidenceGraph(String oneLineFormat, MotifFactory factory, boolean lineContainsPermKey){

		
		Scanner scan = GeneralToolbox.generateScanner(oneLineFormat);

		if (lineContainsPermKey){
			@SuppressWarnings("unused")
			String permGroup = scan.next();
		}
			
		String s = scan.next();
		
		motif = factory.createMotifFromString(s);
		freqVec = new FreqVec();
		
		//System.out.println(s);
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			int f = scan.nextInt();
			//System.out.println(f);
			freqVec.setFreq(i,f);
			
			
		}
		
		probVec = new ProbabilityVector();
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			double p = scan.nextDouble();
			//System.out.println(p);
			probVec.setProb(i,p);
		}
		
		scan.close();

		
		
	}
	
	public BLSConfidenceGraph(String oneLineFormat, MotifFactory factory){
		
		Scanner scan = GeneralToolbox.generateScanner(oneLineFormat);
		
		
		@SuppressWarnings("unused")
		String permGroup = scan.next();
		
		String s = scan.next();
		
		motif = factory.createMotifFromString(s);
		freqVec = new FreqVec();
		
		//System.out.println(s);
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			int f = scan.nextInt();
			//System.out.println(f);
			freqVec.setFreq(i,f);
			
			
		}
		
		probVec = new ProbabilityVector();
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			double p = scan.nextDouble();
			//System.out.println(p);
			probVec.setProb(i,p);
		}
		
		scan.close();
	
	}
	
	
	
	//SETTERS

	public void setMotif(Motif motif) {
		this.motif = motif;
	}
	
	public void setFreqVec(FreqVec freqVec) {
		this.freqVec = freqVec;
	}
	
	public void setProbVec(ProbabilityVector probVec) {
		this.probVec = probVec;
	}
	
		
	//GETTERS

	@Override
	public Motif getMotif() {
		return motif;
	}

	@Override
	public int getFreq(int i){
		return freqVec.getFreq(i);
	}

	@Override
	public double getProbValue(int i) {
		return probVec.get(i).getValue();
	}
	
	@Override
	public FreqVec getFreqVec() {
		return freqVec;
	}
	
	//METHODS
	
	@Override
	public String toString(){
		StringBuilder output = new StringBuilder();
		output.append("Motif="+motif.toString());
		output.append(newline);
		int [] thresholds=BLS.getBLSThresholds();
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			output.append(thresholds[i]);
			output.append(tab);
			output.append(probVec.get(i));
			output.append(tab);
			output.append(freqVec.getFreq(i));
			output.append(tab);
			output.append(newline);


		}
		return output.toString();
	}
	
	public String toOneLineStringFormat(){
		StringBuilder output = new StringBuilder();
		output.append(motif.toString());
		output.append(tab);

		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			output.append(freqVec.getFreq(i));
			output.append(tab);
		}
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			output.append(probVec.get(i));
			output.append(tab);
		}
		return output.toString();
	}
	
	@Override
	public void toText(Text t){
		t.set(this.toOneLineStringFormat());
	}
	
	
	





	


	
	
	

	
}
