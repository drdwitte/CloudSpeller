package be.iminds.cloudspeller.driver;

import be.iminds.cloudspeller.motifalgorithms.DiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.DiscoveryAlgorithmFactory;
import be.iminds.cloudspeller.motifalgorithms.MotifAlgType;
import be.iminds.cloudspeller.motifalgorithms.MotifSearchSpace;
import be.iminds.cloudspeller.motifalgorithms.BenchmarkType;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.MMFactory;
import be.iminds.cloudspeller.motifmodels.MotifFactory;
import be.iminds.cloudspeller.motifmodels.MotifFreqVec;

import be.iminds.cloudspeller.motifpermutationgroups.IUPACContentFactory;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContentFactory;
import be.iminds.cloudspeller.motifpermutationgroups.MotifPermutationGroup;

import org.apache.commons.lang.NotImplementedException;
import org.apache.hadoop.conf.Configuration;

import be.iminds.cloudspeller.output.BLSConfGraphFactory;
import be.iminds.cloudspeller.output.ConfidenceGraphFactory;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.SimultaneousOccurrenceAndConfidenceFiltering;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculatorFactory;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculatorFactory;

import be.iminds.cloudspeller.alphabets.BasePairAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.IndexStructureFactory;
import be.iminds.cloudspeller.indexing.NodeDecorationFactory;
import be.iminds.cloudspeller.indexing.SeqIDDecorationFactory;
import be.iminds.cloudspeller.indexing.BitSetDecorationFactory;
import be.iminds.cloudspeller.input.GeneFamily;


public class CloudSpeller {

	private Configuration conf;
	
	private final String [] parameters  = {"MotifAlgorithm_Type", "Index_Structure", 
			"Node_Decoration_Type", "ConservationScore", "BLS_Thresholds", "Kmin", "Kmax", "Max_Degenerate_Positions", 
			"Motif_Alphabet" , "Filter_Type", 
			"Family_Cutoff", "Confidence_Cutoff", "Background_Group_Size"};
	
	//mapper
	private MotifAlgType motifAlgType;
	private IndexStructureFactory indexStructureFactory;
	private ConservationScoreCalculatorFactory calculatorFactory;
	private MotifSearchSpace motifSearchSpace;
	private MotifFactory motifFactory;

	//reducer
	private MotifContentFactory motifContentFactory;

	//be.iminds.cloudspeller.output
	private ConfidenceGraphFactory confGraphFactory;
	private ConfidenceGraphRestrictions confidenceGraphRestrictions;

	boolean mapperInitialized=false;
	boolean reducerInitialized=false;
	boolean combinerInitialized=false;
	
	/**
	 * Initialize all factories based on Settings in JobConf
	 * 
	 */
	public CloudSpeller(Configuration conf) {
		this.conf=conf;
	}

	
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
	
		for (String parameter : parameters){
			sb.append(parameter+"="+conf.get(parameter)+"\n");
		}
		

		
		return sb.toString();
		
	}
	
	/**
	 * Lazy initialization (only when required)
	 * Choose motif algorithm, index structure, conservation score, number of threshold
	 * intervals, set motif search space
	 * @param conf
	 */
	public void initializeMapperFactories(){
		
		if (mapperInitialized){
			return;
		}
		initializeMotifAlgorithmType();
		initializeMotifContentType(); 
		initializeIndexStructure();
		initializeConservationScore();
		initializeMotifSearchSpace();
		mapperInitialized=true;
	}	
	

	
	public void initializeCombinerFactories() {
		
		if (combinerInitialized)
			return;
		
		initializeMotifAlgorithmType(); //to deserialize MotifFreq
		initializeMotifContentType(); //to deserialize Key
		initializeConservationScore(); //
		
		combinerInitialized=true;
	}
	
	public void initializeReducerFactories(){
		if (reducerInitialized)
			return;
		
		initializeMotifAlgorithmType(); //to deserialize MotifFreq
		initializeMotifContentType(); //to deserialize Key
		initializeConservationScore(); //
		initializeOutputType();
	
		reducerInitialized=true;
		
	}
	
	private void initializeMotifAlgorithmType(){
		
		//Motif algorithm type: EXACT/INEXACT/FAKE + set motif factory 
		
		motifAlgType=null;
		String temp=conf.get("MotifAlgorithm_Type","EXACT");
		if (temp.contains("EXACT")) {
			String strAlg = temp;
			
			temp =conf.get("Motif_Alphabet", "TWOFOLDSANDN");
			IUPACType type;
			if (temp.equals("BASEPAIRS")) {
				type = IUPACType.BASEPAIRS;
			} else if (temp.equals("DONTCARES")) {
				type = IUPACType.DONTCARES;
			} else if (temp.equals("TWOFOLDSANDN")) {
				type = IUPACType.TWOFOLDSANDN;
			} else if (temp.equals("FULL")) {
				type = IUPACType.FULL;
			} else {
				throw new NotImplementedException("IUPAC type not supported");
			}
			motifFactory = new IUPACFactory(type);
			
			if (strAlg.equals("EXACT")){
			
				motifAlgType = MotifAlgType.EXACT;
				
			} else if (strAlg.equals("EXACTAB")){
			
				motifAlgType= MotifAlgType.EXACTAB;
			} else {
				throw new NotImplementedException(strAlg+" not supported: choose EXACT or EXACTAB");
			}
			
					
					 
			
			
			
			
		} else if (temp.equals("INEXACT")) {
			motifAlgType = MotifAlgType.INEXACT;
			motifFactory = new MMFactory();

		} else {
			motifAlgType = MotifAlgType.FAKE;
			
			temp =conf.get("Motif_Alphabet", "TWOFOLDSANDN");
			IUPACType type;
			if (temp.equals("BASEPAIRS")) {
				type = IUPACType.BASEPAIRS;
			} else if (temp.equals("DONTCARES")) {
				type = IUPACType.DONTCARES;
			} else if (temp.equals("TWOFOLDSANDN")) {
				type = IUPACType.TWOFOLDSANDN;
			} else if (temp.equals("FULL")) {
				type = IUPACType.FULL;
			} else {
				throw new NotImplementedException("IUPAC type not supported");
			}
			motifFactory = new IUPACFactory(type);
			
			String bType = conf.get("Benchmark_Type","EMITRANDOM");
			BenchmarkType benchmarkType;
			if (bType.equals("EMITSINGLEMOTIF")){
				System.out.println("BM: Emit single motif");
				benchmarkType = BenchmarkType.EMITSINGLEMOTIF;

			} else if (bType.equals("EMITSINGLEPERMGROUP")){
				System.out.println("BM: Emit single permutation group");
				benchmarkType = BenchmarkType.EMITSINGLEPERMGROUP;
			} else {
				System.out.println("BM: Emit random motifs");
				benchmarkType = BenchmarkType.EMITRANDOM;
			}
			
			int maximalNumberOfMotifs = conf.getInt("Maximal_Number_Of_Motifs",1000);
			System.out.println("Number of motifs tried: " +maximalNumberOfMotifs);
			DiscoveryAlgorithmFactory.setBenchmarkType(benchmarkType);
			DiscoveryAlgorithmFactory.setMaximalNumberOfMotifs(maximalNumberOfMotifs);
		}
		
		int maxMotifLength = conf.getInt("Kmax",6);
		motifFactory.setMaxLength(maxMotifLength); //required for compression
		MotifFreqVec.setMotifFactory(motifFactory);
	}
	
	private void initializeMotifContentType() {
		if (motifAlgType==null){
			initializeMotifAlgorithmType();
		}
		
		motifContentFactory=null;
		switch(motifAlgType){
		
		case FAKE:
			motifContentFactory = new IUPACContentFactory(motifFactory);
			break;
		case EXACT:
			motifContentFactory = new IUPACContentFactory(motifFactory);
			break;
		case EXACTAB:
			motifContentFactory = new IUPACContentFactory(motifFactory);
			break;
		case INEXACT:
			throw new NotImplementedException();
			//break;
		default: 
			throw new NotImplementedException();
		}
		
		
	}
	
	private void initializeIndexStructure(){

		indexStructureFactory=null;
		String temp=conf.get("Index_Structure","GST");
		if (temp.equals("SA")){
			//indexStructureFactory = new SAFactory();
			throw new NotImplementedException();
		} else if (temp.equals("GST")){
			boolean withRevComp=conf.getBoolean("Reverse_Complements",true);
			int maxTreeDepth = conf.getInt("Kmax",6);
			temp=conf.get("Node_Decoration_Type","BITS");
			NodeDecorationFactory nodeDecoFac;
			if (temp.equals("SETS")){
				nodeDecoFac = new SeqIDDecorationFactory();
			} else if (temp.equals("BITS")) {
				nodeDecoFac = new BitSetDecorationFactory();
			} else {
				throw new NotImplementedException();
			}
			indexStructureFactory = new GSTFactory(maxTreeDepth,withRevComp,nodeDecoFac);
			
		
		} else {
			//indexStructureFactory = new FISFactory();
		}
	}
	
	private void initializeConservationScore(){
		//Choose conservation score calculator and set number of frequency intervals
		String temp;
		calculatorFactory=null;
		temp=conf.get("ConservationScore","BLS");
		
					
		if (temp.equals("BLS")){
			String thresholdString = conf.get("BLS_Thresholds","");
			
			if (thresholdString.length()>0){ //check if set
				String [] strThresholds = thresholdString.split(",");
				
				int [] thresholds = new int[strThresholds.length];
				for (int i=0; i<thresholds.length; i++){
					thresholds[i] = Integer.parseInt(strThresholds[i]);
				}
				
				//check if sorted
				boolean badArray=false;
				for (int i=1; i<thresholds.length; i++){
					if (thresholds[i]<=thresholds[i-1]){
						badArray=true;
						break;
					}
				}
				
				if (!badArray){
					BLS.initializeBLSConstants(thresholds);
				} else {
					System.err.println("Bad BLS array using default [10,30,50,70,90]");
				}
			}
			
			
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			calculatorFactory = new BLSCalculatorFactory(new BLS(BLS.MIN));
		} else if (temp.equals("SC")){
			throw new NotImplementedException();
		} else if (temp.equals("Q")){ 
			throw new NotImplementedException();
		} else { 
			throw new NotImplementedException();
		}
	}



	private void initializeOutputType()  {
		
		String filterType=conf.get("Filter_Type","SIMULTANEOUS");
		int famCutoff=conf.getInt("Family_Cutoff",5);
		int confCutoff=conf.getInt("Confidence_Cutoff",90);
		int permutationGroupSize=conf.getInt("Background_Group_Size",1000);
		
		MotifPermutationGroup.setBackgroundGroupSize(permutationGroupSize);
				
		if (filterType.equals("SIMULTANEOUS")){
			confidenceGraphRestrictions = 
				new SimultaneousOccurrenceAndConfidenceFiltering(confCutoff,famCutoff);
		} else {
			throw new NotImplementedException();
		}
		
		String temp=conf.get("ConservationScore","BLS");
		if (temp.equals("BLS")){
			confGraphFactory = new BLSConfGraphFactory();
		} else {
			throw new NotImplementedException();
		}

		
	}
	
	private void initializeMotifSearchSpace(){
		//Choose motif search space (kmin, kmax, maxDegeneracy, alphabet)
		motifSearchSpace=null;
		int kmin=conf.getInt("Kmin",6);
		int kmax=conf.getInt("Kmax",6);
		int maxDegPositions=conf.getInt("Max_Degenerate_Positions",0);
		//System.out.println("searchSpace init: "+kmin+" "+kmax+" "+maxDegPositions);
		String temp=conf.get("Motif_Alphabet","BASEPAIRS");
		if (temp.equals("BASEPAIRS")){
			motifSearchSpace = new MotifSearchSpace(kmin,kmax,0,new BasePairAlphabet());
		} else if (temp.equals("DONTCARES")){
			motifSearchSpace = new MotifSearchSpace(kmin,kmax,maxDegPositions,new IUPACAlphabet(IUPACType.DONTCARES));
		} else if (temp.equals("TWOFOLDSANDN")){
			motifSearchSpace = new MotifSearchSpace(kmin,kmax,maxDegPositions,new IUPACAlphabet(IUPACType.TWOFOLDSANDN));
		} else if (temp.equals("IUPAC")){
			motifSearchSpace = new MotifSearchSpace(kmin,kmax,maxDegPositions,new IUPACAlphabet(IUPACType.FULL));
		} else {
			throw new NotImplementedException();
		}
		
	}
	
	//GETTERS
	
	public DiscoveryAlgorithm getDiscoveryAlgorithm(GeneFamily geneFamily) {
		DiscoveryAlgorithm alg = DiscoveryAlgorithmFactory.createAlgorithm(geneFamily.getSequences(),motifAlgType);
		alg.setDataStructure(indexStructureFactory);
		alg.setConservationScoreCalculator(calculatorFactory.createCalculator(geneFamily));
		alg.setSearchSpace(motifSearchSpace);
		return alg;
	}

	public ConfidenceGraphFactory getConfGraphFactory() {
		return confGraphFactory;
	}

	public ConfidenceGraphRestrictions getConfGraphRestrictions() {
		return confidenceGraphRestrictions;
	}
	
	public MotifContentFactory getMotifContentFactory() {
		return motifContentFactory;
	}

	public MotifFactory getMotifFactory() {
		return motifFactory;
	}



	public MotifAlgType getMotifAlgType() {
		return motifAlgType;
	}

}

