package be.iminds.cloudspeller.postprocessing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.SortedSet;
import java.util.TreeSet;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;


import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.phylogenetics.BLS;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;


/**
 * Extract motifs with the highest number of targets
 * @author ddewitte
 *
 */

public class BestMotifsExtractor {
	
	private String filename;
	private int numberOfMotifs;
	private final int multiplicFactorForOptimization=10;
	private Alphabet alphabet;
	private IUPACFactory factory; 
	private SortedSet<MotifOcc> bestMotifs = new TreeSet<MotifOcc>();
	private MotifAligner aligner;
	
	
	public BestMotifsExtractor(String filename, int numberOfMotifs, IUPACType type){
		this.filename = filename;
		this.numberOfMotifs = numberOfMotifs;
		this.factory = new IUPACFactory(type);
		this.alphabet = new IUPACAlphabet(type);
		this.aligner= new MotifAligner(this.alphabet, new JacardDistanceCalculator(this.alphabet));
	}
	
	private void printBestMotifs(){ 
		for (MotifOcc mOcc : bestMotifs){
			System.out.println(mOcc);
		}
	}
	
	
	public void filterConfidenceGraphs(ConfidenceGraphRestrictions restrictions, String output) throws IOException{
		
		BufferedReader in = new BufferedReader(new FileReader(filename));
		BufferedWriter out = new BufferedWriter(new FileWriter(output));
		
	
		String oneLine;
		//int [] bls = BLS.getBLSThresholds();
		int numberOfLinesRead=0;
		
		while ((oneLine = in.readLine())!= null){
		
			BLSConfidenceGraph graph = new BLSConfidenceGraph(oneLine,factory);
		
			if (restrictions.checkRestrictions(graph)){
				out.write(graph.toOneLineStringFormat());
			}
			
			if (++numberOfLinesRead%1000==0){
				System.out.println("Number of motifs read: " + numberOfLinesRead);
			}
		
		}//outerwhileloop
		in.close();
		out.close();
		
		
	}
	
	public void extractBestMotifsNoDistanceFilter(double confidenceCutoff) throws IOException {
		
		bestMotifs.clear();
		BufferedReader in = new BufferedReader(new FileReader(filename));
		String oneLine;
		int [] bls = BLS.getBLSThresholds();
		int numberOfLinesRead=0;
		while ((oneLine = in.readLine())!= null){
		
			BLSConfidenceGraph graph = new BLSConfidenceGraph(oneLine,factory);
		
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				if (graph.getProbValue(i)>=confidenceCutoff){
					String mStr = graph.getMotif().toString();
					tryToInsertMotif(mStr, bls[i], graph.getFreq(i));
					break;
				}
			
			}//innerforloop
			
			if (++numberOfLinesRead%1000==0){
				System.out.println("Number of motifs read: " + numberOfLinesRead);
			}
		
		}//outerwhileloop
		
		in.close();
		
		SortedSet<MotifOcc> bestMotifsRestr = new TreeSet<MotifOcc>();
		
		for (MotifOcc mOcc : bestMotifs){
			if (mOcc.getMotif().startsWith("N"))
				continue;
			if (mOcc.getMotif().endsWith("N"))
				continue;
			
			bestMotifsRestr.add(mOcc);
				
					
				
			if (bestMotifsRestr.size()==numberOfMotifs){
					break;
			}
			
		}
		
		bestMotifs = bestMotifsRestr;
	}
	
	public void extractBestMotifs(double confidenceCutoff, double minDist) throws IOException{
		
		
		bestMotifs.clear();
		BufferedReader in = new BufferedReader(new FileReader(filename));
		String oneLine;
		int [] bls = BLS.getBLSThresholds();
		int numberOfLinesRead=0;
		while ((oneLine = in.readLine())!= null){
		
			BLSConfidenceGraph graph = new BLSConfidenceGraph(oneLine,factory);
		
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				if (graph.getProbValue(i)>=confidenceCutoff){
					String mStr = graph.getMotif().toString();
					tryToInsertMotif(mStr, bls[i], graph.getFreq(i));
					break;
				}
			
			}//innerforloop
			
			if (++numberOfLinesRead%1000==0){
				System.out.println("Number of motifs read: " + numberOfLinesRead);
			}
		
		}//outerwhileloop
		
		in.close();
		
		SortedSet<MotifOcc> bestMotifsAtMinDistance = new TreeSet<MotifOcc>();
		
		for (MotifOcc mOcc : bestMotifs){
			if (mOcc.getMotif().startsWith("N"))
				continue;
			if (mOcc.getMotif().endsWith("N"))
				continue;
			
			double minDistForThisMotif = calcMinDistance(mOcc, bestMotifsAtMinDistance);
			
			if (minDistForThisMotif >= minDist ){
				bestMotifsAtMinDistance.add(mOcc);
				
					
				
				if (bestMotifsAtMinDistance.size()==numberOfMotifs){
					break;
				}
			}
		}
		
		bestMotifs = bestMotifsAtMinDistance;
		
	}
	
	public void printToFile(String filename) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		for (MotifOcc mOcc : bestMotifs){
			out.write(mOcc.toString());
			out.write("\n");
		}
		out.close();
	}

	private double calcMinDistance(MotifOcc currentOccurrence, SortedSet<MotifOcc> bestMotifsAtMinDistance) {
		
		double minDistance=Double.MAX_VALUE;
		
		for (MotifOcc mOcc : bestMotifsAtMinDistance){
			double newDist = aligner.calculateSimilarityDistance(mOcc.getMotif()
													,currentOccurrence.getMotif());
			
			if (newDist < minDistance)
				minDistance = newDist;
		}
		return minDistance;
		
	}

	private void tryToInsertMotif(String mStr, int bls, int freq) {
		if (bestMotifs.size()==multiplicFactorForOptimization*numberOfMotifs){
			if (bestMotifs.last().getFamOcc() < freq){
				bestMotifs.remove(bestMotifs.last());
				bestMotifs.add(new MotifOcc(mStr,bls,freq));
			}
		} else {
			bestMotifs.add(new MotifOcc(mStr,bls,freq));
		}
		
		
	}
	
	public static void main (String [] args) throws IOException {
		
		String filename = "outputDCAll.txt";
		int numberOfMotifs=100;
		
		FreqVec.setNumberOfIntervals(6);
		BLS.initializeBLSConstants(40,10,6);
		
		BestMotifsExtractor extractor = new BestMotifsExtractor(filename, numberOfMotifs, IUPACType.DONTCARES);
		
		double minDist = 2.5;
		double confCutoff = 99.0;
		extractor.extractBestMotifs(confCutoff,minDist);
		
		
		extractor.printBestMotifs();
	}

	
}


class MotifOcc implements Comparable<MotifOcc> {
	
	private String motif;
	private double bls;
	private int famOcc;
	
	
	
	public String getMotif() {
		return motif;
	}

	public int getFamOcc() {
		return famOcc;
	}
	
	public double getBLS(){
		return bls;
	} 
	
	public String toString(){
		return motif+"\t"+bls+"\t"+famOcc;
	}

	
	
	/**
	 * Constructor
	 * @param m
	 * @param bls
	 * @param famOcc
	 */
	public MotifOcc(String m, int bls, int famOcc){
		this.motif = m;
		this.bls = bls;
		this.famOcc=famOcc;
	}


	@Override
	/**
	 * if m1 has higher famOcc than m2 then m1<m2 => comparison is occ2-occ1
	 */
	public int compareTo(MotifOcc o) {
		return o.famOcc - this.famOcc;
	}
}





	
	
	
	
		
	
	
	
