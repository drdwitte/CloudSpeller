package be.iminds.cloudspeller.motifpermutationgroups;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.driver.OutputExtractor;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFreqVec;

import be.iminds.cloudspeller.output.ConfidenceGraph;
import be.iminds.cloudspeller.output.ConfidenceGraphFactory;
import be.iminds.cloudspeller.output.ProbabilityVector;




public class MotifPermutationGroup {

	private static int backgroundGroupSize=1000;
	private MotifContent content;
	private double [] backgroundOccs;
	private Map<Motif,FreqVec> motifFamOccMap= new HashMap<Motif,FreqVec>();
	
	/**
	 * Constructor
	 * @param content
	 */
	public MotifPermutationGroup(MotifContent content){
		this.content=content;
	}

	
	//SETTERS
	
	public static void setBackgroundGroupSize(int size){
		backgroundGroupSize=size;
	}
	
	//GETTERS
	
	public double [] getBackgroundModel(){
		return backgroundOccs;
	}
	
	public Map<Motif,FreqVec> getAggregatedMotifMap(){
		return motifFamOccMap;
	}
		
	//METHODS
	
	public void addMotifFreq(MotifFreqVec mFreq){
		FreqVec f=motifFamOccMap.get(mFreq.getMotif());
		if (f!=null){
			f.add(mFreq.getVec());
		} else {
			motifFamOccMap.put(mFreq.getMotif(),mFreq.getVec());
		}
	}
	
	public void generateBackgroundModel(){
		
		
		Set<Motif> backgrSet=content.createPermutationGroup(backgroundGroupSize);
		
		backgroundOccs = new double [FreqVec.getNumberOfIntervals()];
		final int noMotifOccurrence=0;
		
		ArrayList< ArrayList<Integer> > occs= 
			new ArrayList< ArrayList<Integer> >(FreqVec.getNumberOfIntervals());

		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			occs.add(new ArrayList<Integer>());
		}
		
		Iterator<Motif> iterator=backgrSet.iterator();
		while (iterator.hasNext()){
			Motif motif=iterator.next();
			FreqVec fV=motifFamOccMap.get(motif);
			if (fV!=null){
				for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
					occs.get(i).add(fV.getFreq(i));
					
				}
			} else {
				for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
					occs.get(i).add(noMotifOccurrence);
				}
			}
		}
		
		
		for (int i=0; i<backgroundOccs.length; i++){
			backgroundOccs[i]=findMedian(occs.get(i));
		}
	
	}
		
	public static double findMedian(ArrayList<Integer> a){
		
		Collections.sort(a);
		int numEl=a.size();
		if (numEl<2){
			if (numEl==0){
				return 0;
			} else {
				return a.get(0);
			}
		}
		else if (numEl%2==0){
			int i=numEl/2-1;
			return 0.5*(a.get(i)+a.get(i+1));
		} else {
			int i=numEl/2;
			return a.get(i);
		}
	}
	
	public void generateConfidenceGraphs(ConfidenceGraphFactory confGraphFactory, OutputExtractor extractor){
		
		for (Map.Entry<Motif,FreqVec> mFreqVec : motifFamOccMap.entrySet()){
			Motif m= mFreqVec.getKey();
			FreqVec f=mFreqVec.getValue();
			ConfidenceGraph graph=confGraphFactory.createGraph(m,f,new ProbabilityVector(f, backgroundOccs));
			extractor.extract(graph);
		}
	}

	public int groupSize() {
		return motifFamOccMap.size();
	}
	
	
	public double [] calculateMean(Set<Motif> backgrSet){
		double  [] mean = new double [FreqVec.getNumberOfIntervals()];
		int N = backgrSet.size();
		
		for (Motif m : backgrSet){
			
			FreqVec fV=motifFamOccMap.get(m);
			if (fV!=null){
				for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
					mean[i]+=fV.getFreq(i);
					
				}
			}
		
		}
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			mean[i]/=N;
			
		}
		
		return mean;
	}
	
	public double [] calculateMeanSquare(Set<Motif> backgrSet){
		double  [] mean = new double [FreqVec.getNumberOfIntervals()];
		int N = backgrSet.size();
		
		for (Motif m : backgrSet){
			
			FreqVec fV=motifFamOccMap.get(m);
			if (fV!=null){
				for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
					mean[i]+=fV.getFreq(i)*fV.getFreq(i);
					
				}
			}
		
		}
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			mean[i]/=N;
			
		}
		
		return mean;
	}
	
	
	/**
	 * Warning underflow risk! (probably not)
	 * Square root of variance
	 * @param mean
	 * @param meanSquare
	 * @return
	 */
	public double [] calculateVariation(double [] mean, double [] meanSquare){
		
		double  [] var = new double [FreqVec.getNumberOfIntervals()];
		
		for (int i=0; i<var.length; i++){
			var[i]=Math.sqrt(meanSquare[i]-mean[i]*mean[i]);
		}
		
		return var;
	}
	
	public int [] calculatePValueFamOcc(double[] alpha, int freqIndex, Set<Motif> backgrSet){
		int  [] pvalueOcc = new int [alpha.length];
		
		List<Integer> occList = new ArrayList<Integer>();
		
		
		for (Motif m : backgrSet){
			
			FreqVec fV=motifFamOccMap.get(m);
			if (fV!=null){
				occList.add(fV.getFreq(freqIndex));
			} else {
				occList.add(0);
			}
		
		}
		
		Collections.sort(occList);
		
		for (int i=0; i<alpha.length; i++){
			double alph = alpha[i];
			
			int elsRight = (int)Math.ceil(alph * occList.size());
			int index = occList.size()-elsRight;
			
			
			//shift to F value which is enough to guarantee pvalue
			//example 1111 22 33 4 55 66 7 pvalue 2/N can only be obtained for F=7! since 
			//F=6 is for 3/N!!
			
			if (index>0){
				while (occList.get(index) == occList.get(index-1) && index < occList.size()-1){
					index++;
				}
			}
			pvalueOcc[i]=occList.get(index);
						
		}
		
		
		return pvalueOcc;
	}


	public double [] getPercentageInSigmaInterv(int numberOfSigma, Set<Motif> backgrSet, double[] mean, double[] sigma) {
		
		double [] percentages = new double[mean.length];
		
		for (int i=0; i<percentages.length; i++){
			percentages[i] = getPercInSigma(numberOfSigma, backgrSet, mean[i], sigma[i],i);
		}
		
		return percentages;
		
		
	}
	
	private double getPercInSigma(int nSigma, Set<Motif> backgrSet, double mean, double sigma, int freqIndex) {
		
		int N = backgrSet.size();
		int numberInBulk=0;
		for (Motif m : backgrSet){
			FreqVec fV=motifFamOccMap.get(m);
			int freq = 0;
			if (fV!=null){
				freq = fV.getFreq(freqIndex);
			} 
			
			if (freq <= mean+nSigma*sigma  && freq >= mean-nSigma*sigma){
				numberInBulk++;
			} 
		}
		
		double percentage = (1.0*numberInBulk)/N;
		return percentage;
		
		
	}
	
}
