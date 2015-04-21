package be.iminds.cloudspeller.input;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifpermutationgroups.IUPACContent;


public class ShuffledDatasetPreparator extends DatasetPreparator {

	
	public static void main (String [] args ) throws IOException {
		
		//runShuffle(1);
		runFullRandomShuffle();
	}
	
	public static void runFullRandomShuffle() throws IOException{
		
		Map<Gene,String> allGeneSeqs=readAllSequenceFiles("");
		Map<String,Set<Gene>> geneFams = readSelectionFile("iORTHO_2evid_osa.selection",genePrefixOrgMap);
		
		System.out.println("Number of gene families in selection file: "+geneFams.size());
		
		GeneFamily.setGeneralNewick(readNewickFile("speciesTree.txt"));
			
		int cnt=0;
		for (Map.Entry<Gene,String> e : allGeneSeqs.entrySet()){
			
			cnt++;
			if (cnt%1000==0){
				System.out.println(cnt+"/"+allGeneSeqs.size());
			}
			
			List<Character> chars = IUPACContent.charListFromString(e.getValue());
			
			Collections.shuffle(chars);
			
			StringBuilder newSeq = new StringBuilder();
			
			for (Character c : chars){
				newSeq.append(c);
			}
			
			e.setValue(newSeq.toString());
			
			newSeq=null;
		}
		
		
		
		
		Set<GeneFamily> geneFamilies = createGeneFamilies(geneFams, allGeneSeqs);

		System.out.println("Printing to file...");
		int numFamPerFile=10;
		
		String path = "/home/ddewitte/Desktop/CloudSpellerExperimentsFinal/RandomDatasets/";
		String map = "MonocotFamiliesAgg"+numFamPerFile+"PBL1em6Shuffled1/";
		new File(path+map).mkdir();
			
		printGeneFamiliesInGroups(path+map,geneFamilies,numFamPerFile);
		System.out.println("Done! Bye!");
	}
	
	public static void runShuffle(int shuffleLength) throws IOException{
		int seqLength = 2000;
		int numFrags = seqLength/shuffleLength;
		
		Map<Gene,String> allGeneSeqs=readAllSequenceFiles("");
		Map<String,Set<Gene>> geneFams = readSelectionFile("iORTHO_2evid_osa.selection",genePrefixOrgMap);
		
		System.out.println("Number of gene families in selection file: "+geneFams.size());
		System.out.println("Shuffle length: "+shuffleLength);
		
		GeneFamily.setGeneralNewick(readNewickFile("speciesTree.txt"));
			
		Map<Gene,String> zmaGeneSequencesInDataset = extractAllGenesWithPrefixInDataset("ZM",allGeneSeqs,geneFams);
		Map<Gene,String> bdiGeneSequencesInDataset = extractAllGenesWithPrefixInDataset("BD",allGeneSeqs,geneFams);
		Map<Gene,String> sbiGeneSequencesInDataset = extractAllGenesWithPrefixInDataset("SB",allGeneSeqs,geneFams);
		Map<Gene,String> osaGeneSequencesInDataset = extractAllGenesWithPrefixInDataset("OS",allGeneSeqs,geneFams);
		
		System.out.println("zmaOrig");
		analyzeKmerDistribution(zmaGeneSequencesInDataset);
		System.out.println("bdiOrig");
		analyzeKmerDistribution(bdiGeneSequencesInDataset);
		System.out.println("sbiOrig");
		analyzeKmerDistribution(sbiGeneSequencesInDataset);
		System.out.println("osaOrig");
		analyzeKmerDistribution(osaGeneSequencesInDataset);
				
		reshuffle(zmaGeneSequencesInDataset, numFrags);
		reshuffle(bdiGeneSequencesInDataset, numFrags);
		reshuffle(sbiGeneSequencesInDataset, numFrags);
		reshuffle(osaGeneSequencesInDataset, numFrags);
		
		System.out.println("zmaShuff");
		analyzeKmerDistribution(zmaGeneSequencesInDataset);
		System.out.println("bdiShuff");
		analyzeKmerDistribution(bdiGeneSequencesInDataset);
		System.out.println("sbiShuff");
		analyzeKmerDistribution(sbiGeneSequencesInDataset);
		System.out.println("osaShuff");
		analyzeKmerDistribution(osaGeneSequencesInDataset);
		
		updateGeneSequences(allGeneSeqs, zmaGeneSequencesInDataset);
		updateGeneSequences(allGeneSeqs, bdiGeneSequencesInDataset);
		updateGeneSequences(allGeneSeqs, sbiGeneSequencesInDataset);
		updateGeneSequences(allGeneSeqs, osaGeneSequencesInDataset);
		
		
		
		Set<GeneFamily> geneFamilies = createGeneFamilies(geneFams, allGeneSeqs);
			
		
		System.out.println("Printing to file...");
		int numFamPerFile=10;
		
		String path = "/home/ddewitte/Bureaublad/CloudSpellerExperiments/RandomDatasets/";
		String map = "MonocotFamiliesAgg"+numFamPerFile+"PBL1em6Shuffled"+shuffleLength+"/";
		new File(path+map).mkdir();
			
		printGeneFamiliesInGroups(path+map,geneFamilies,numFamPerFile);
		System.out.println("Done! Bye!");
	}
	


	private static Map<Gene, String> extractAllGenesWithPrefixInDataset(
			String prefix, Map<Gene, String> allGeneSeqs, Map<String, Set<Gene>> geneFams) {
		
		
		Map<Gene,String> genesWithPrefix = new HashMap<Gene, String>();
		
		Set<Gene> genes = getAllGenesInDataset(prefix,geneFams);
		
		for (Gene g : genes){
			String seq = allGeneSeqs.get(g);
			if (seq!=null){
				genesWithPrefix.put(g,seq);
			}
		}
		return genesWithPrefix;
	}

	private static Set<Gene> getAllGenesInDataset(String prefix, Map<String, Set<Gene>> geneFams) {
		
		Set<Gene> genes = new HashSet<Gene>();
		for (Map.Entry<String,Set<Gene>> e : geneFams.entrySet()){
			for (Gene g : e.getValue()){
				if (g.getOrganism().equals(prefix)){
					genes.add(g);
				}
			}
		}
		
		return genes;
		
	}

	private static void updateGeneSequences(Map<Gene, String> allGeneSeqs,
			Map<Gene, String> newGeneSequences) {
		
		allGeneSeqs.putAll(newGeneSequences);
		
		
	}


	private static void reshuffle(
			Map<Gene, String> geneSequencesInDataset, int numFrags) {
		
		
		List<String> fragmentList = createFragments(geneSequencesInDataset,numFrags);
		Collections.shuffle(fragmentList);
		reassembleGeneSeqMap(geneSequencesInDataset,fragmentList);
		
		
		
	}

	private static void reassembleGeneSeqMap(
			Map<Gene, String> geneSequencesInDataset, List<String> fragmentList) {
		
		int fragsPerSeq = fragmentList.size()/geneSequencesInDataset.size();
		int fragCounter = 0;
		
		for (Map.Entry<Gene,String> e : geneSequencesInDataset.entrySet()){
			StringBuilder newSeq = new StringBuilder();
			
			for (int i=0; i<fragsPerSeq; i++){
				newSeq.append(fragmentList.get(fragCounter++));
			}
			e.setValue(newSeq.toString());
		}
		
		
		
		
		
		
	}

	/**
	 * Currently delta not used
	 * 
	 */
	private static List<String> createFragments(
			Map<Gene, String> geneSequencesInDataset, int numFrags) {
			
		List<String> frags = new ArrayList<String>();
		for (Map.Entry<Gene,String> e : geneSequencesInDataset.entrySet()){
			
			String seq = e.getValue();
			
			int fragLength = seq.length()/numFrags;
			String frag;
			for (int i=0; i<numFrags-1; i++){
				frag = seq.substring(i*fragLength,(i+1)*fragLength);
				frags.add(frag);
			}
			
			frag = seq.substring((numFrags-1)*fragLength);
			frags.add(frag);
		}
		
		return frags;
		
	}

	/**
	 * Analyze 3mer distribution
	 * @param zmaGeneSequencesInDataset
	 */
	private static void analyzeKmerDistribution(
			Map<Gene, String> geneSequencesInDataset) {
		
		String alph = "ACGT";
		int kmerLength = 3;
		
		Map<String,Integer> distribution = new HashMap<String,Integer>();
		
		//intialize
		for (int i=0; i<alph.length(); i++){
			for (int j = 0; j < alph.length(); j++) {
				for (int k = 0; k < alph.length(); k++) {
					distribution.put(""+alph.charAt(i)+alph.charAt(j)+alph.charAt(k),0);
				}
			}
		}
		int counter=0;
		for (Map.Entry<Gene,String> e : geneSequencesInDataset.entrySet()){
			counter++;
			if (counter%1000==0){
				System.out.print("*");
			}
			
			String seq = e.getValue();
			for (int i=0; i<(seq.length()+1-kmerLength); i++){
				String kmer = seq.substring(i,i+kmerLength);
				String revKmer = (new IUPACMotif(kmer,0)).getComplement().toString();
				
				
				
				Integer freq1 = distribution.get(kmer);
				
				if (freq1==null){ //contains Ns which are masked in discovery
					continue;
				}
				
				distribution.put(kmer,freq1+1);
				distribution.put(revKmer,freq1+1);

				
			}
			
		}
		
		
		System.out.println();
		//print
		for (int i=0; i<alph.length(); i++){
			for (int j = 0; j < alph.length(); j++) {
				for (int k = 0; k < alph.length(); k++) {
					String key = ""+alph.charAt(i)+alph.charAt(j)+alph.charAt(k);
					int freq = distribution.get(key);
					System.out.println(key + "\t" + freq);
				}
			}
		}
		
		
		
		
	}
}
