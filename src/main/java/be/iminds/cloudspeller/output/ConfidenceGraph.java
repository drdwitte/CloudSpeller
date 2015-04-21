package be.iminds.cloudspeller.output;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;

import org.apache.hadoop.io.Text;


public interface ConfidenceGraph {

	//GETTERS
	public Motif getMotif();
	int getFreq(int i);
	public FreqVec getFreqVec();
	double getProbValue(int i);
	
	//METHODS
	public String toString();
	public void toText(Text t);
	
	

}
