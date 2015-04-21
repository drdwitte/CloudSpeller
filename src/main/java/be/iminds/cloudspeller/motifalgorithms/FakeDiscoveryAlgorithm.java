package be.iminds.cloudspeller.motifalgorithms;

import be.iminds.cloudspeller.indexing.IndexStructureFactory;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.indexing.SequenceIDSet;
import be.iminds.cloudspeller.input.Sequence;

import java.util.ArrayList;
import java.util.Set;

import be.iminds.cloudspeller.driver.MotifExtractor;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;
import be.iminds.cloudspeller.motifpermutationgroups.IUPACContent;


import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculator;


public class FakeDiscoveryAlgorithm implements DiscoveryAlgorithm {
	
	private int numberOfFakeMotifs;
	private ArrayList<Sequence> sequences;
	private MotifSearchSpace motifSearchSpace;
	private ConservationScoreCalculator calculator;
	private MotifExtractor extractor;
	private BenchmarkType benchmarkType;
	/**
	 * Constructor
	 * @param arrayList
	 * @param numberOfFakeMotifs
	 * @param benchmarkType 
	 */
	public FakeDiscoveryAlgorithm(ArrayList<Sequence> arrayList, int numberOfFakeMotifs, BenchmarkType benchmarkType) {
		this.sequences=arrayList;
		this.numberOfFakeMotifs=numberOfFakeMotifs;
		this.benchmarkType=benchmarkType;
	}

	
	//SETTERS
	
	@Override
	public void setDataStructure(IndexStructureFactory indexStructureFactory) {
		
	}
	

	@Override
	public void setSearchSpace(MotifSearchSpace motifSearchSpace) {
		this.motifSearchSpace=motifSearchSpace;
	}


	@Override
	public void setConservationScoreCalculator(ConservationScoreCalculator c) {
		this.calculator=c;
		
	}
	
	//GETTERS

	
	//METHODS
	
	@Override
	public void runDiscovery(MotifFactory motifFactory) {
		switch (benchmarkType){
		case EMITRANDOM:
			generateRandomMotifs(motifFactory);
			break;
		case EMITSINGLEMOTIF:
			generateSingleMotif(motifFactory);
			break;
		case EMITSINGLEPERMGROUP:
			generateSingleGroup(motifFactory);
			break;
		}
		
		
	}
	
	private void generateSingleMotif(MotifFactory motifFactory) {
		//NOTE this motif must have more As than Ts otherwise it is redundant!
		String singleMotif = "ACGTNMRATGCA"; //3 degenerate positions for benchmarking realistic case
		for (int i=0; i<numberOfFakeMotifs; i++) {
			Motif motif = motifFactory.createMotifFromString(singleMotif);
			ConservationScore score = generateRandomScore();
			if (score!=null){
				extractor.add(motif, score);
			}
		}
	}


	private void generateSingleGroup(MotifFactory motifFactory) {
		//NOTE this motif must have more As than Ts otherwise it is redundant!
		String singleMotif = "ACGTNMRATGCA"; //3 degenerate positions for benchmarking realistic case
		Motif motif = motifFactory.createMotifFromString(singleMotif);
		IUPACContent content = new IUPACContent(motif);
		int numMotifsPerShuffleStep=1000;
		int numberOfIterations = numberOfFakeMotifs / numMotifsPerShuffleStep;
		
		for (int i=0; i<numberOfIterations; i++){
			Set<Motif> backgrSet=content.createPermutationGroup(numMotifsPerShuffleStep);
			for (Motif m : backgrSet){
				ConservationScore score = generateRandomScore();
				if (score!=null){
					extractor.add(m, score);
				}
			}
		}
		
		Set<Motif> backgrSetRemainder=content.createPermutationGroup(numberOfFakeMotifs%numMotifsPerShuffleStep);
		for (Motif m : backgrSetRemainder){
			ConservationScore score = generateRandomScore();
			if (score!=null){
				extractor.add(m, score);
			}
		}
	}


	private void generateRandomMotifs(MotifFactory motifFactory){
		for (int i=0; i<numberOfFakeMotifs; i++) {
			
			int randomLength=generateRandomMotifLength();
			IUPACMotif motif = (IUPACMotif)motifFactory.createRandomMotif(randomLength);
			if (motif.numberOfDegPositions()>motifSearchSpace.getMaxNumberOfDegeneratePositions()){
				i--; continue;
			}
			
			ConservationScore score = generateRandomScore();
			if (score!=null){
				extractor.add(motif, score);
			}
		}
	}
	
	private ConservationScore generateRandomScore(){
		NodeDecoration nodeDeco = new SequenceIDSet(sequences.size());
		nodeDeco.setRandomDeco();
		return calculator.calculateScore(nodeDeco);
	}

	private int generateRandomMotifLength(){
		int kmin=motifSearchSpace.getMinLength();
		int kmax=motifSearchSpace.getMaxLength();
		return kmin+(int)(Math.random()*(kmax-kmin+1));
		
	}
	
	
	
	@Override
	public String toString() {
		
		return "I am a fake discovery algorithm";
	}

	@Override
	public void setMotifExtractor(MotifExtractor extractor) {
		this.extractor=extractor;
	}




}
