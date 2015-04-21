package be.iminds.cloudspeller.motifalgorithms;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;



import be.iminds.cloudspeller.alphabets.Alphabet;

import be.iminds.cloudspeller.indexing.GeneralizedSuffixTree;
import be.iminds.cloudspeller.indexing.ISMonkey;
import be.iminds.cloudspeller.indexing.IndexStructure;
import be.iminds.cloudspeller.indexing.IndexStructureFactory;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.BaseSequence;
import be.iminds.cloudspeller.input.Sequence;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

//IMPORTANT reverse motifs and duplicate motifs in strings are taken care of in the MotifExtractor class
//designed for this algorithm: ABInstantEmitter

public class AlBasedExactDiscoveryAlgorithm extends DeNovoExactDiscoveryAlgorithm {

	private IndexStructureFactory iSFac;
	private int prefixLength=3;
	
	public AlBasedExactDiscoveryAlgorithm(ArrayList<Sequence> seqs) {
		super(seqs);
	}
	
	public void setPrefixLength(int prefixLength){
		this.prefixLength=prefixLength;
	}
	
	
	/*for (int i=0; i<sequences.get(0).length(); i++){
        
        List<Suffix> suffixes = new ArrayList<Suffix>();                       
       
        for (int j=0; j<sequences.size(); j++){
                suffixes.add(new Suffix(j,i));
        }
        index = iSFac.createIndexStructureForSuffixes(sequences, suffixes);

        exploreSubtreeFast(index.getExactISMonkey(factory,motifSearchSpace.getMaximumDegeneracy()));
       
}*/
	
	
	@Override
	public void runDiscovery(MotifFactory factory) {
		System.err.println("run ab discovery");
		ArrayList<Sequence> revSequences = BaseSequence.generateReverseComplements(sequences);
		
		//add sentinels and replace non ACGT chars with sentinel
		//not necessary to repeat this for every position + no deep copy of seqs!
		//long start = System.nanoTime();
		ArrayList<Sequence> sequencesPreprocessed = GeneralizedSuffixTree.preprocessSequences(sequences,false);
		ArrayList<Sequence> revSequencesPreprocessed = GeneralizedSuffixTree.preprocessSequences(revSequences,false);
		//System.err.println("Seq preprocessing: "+(System.nanoTime()-start)/1000000);
		
		//start = System.nanoTime();
		ArrayList<IndexStructure> forwardIndexes = generateIndexes(sequencesPreprocessed);
		ArrayList<IndexStructure> reverseIndexes = generateIndexes(revSequencesPreprocessed);
		//System.err.println("Index building: "+(System.nanoTime()-start)/1000000);
		
		
		
		if (extractor!=null){
			
			SortedSet<String> prefixes = generateAllPrefixes();
			//int counter=0;
			//start = System.nanoTime();
			prefixLoop:for (String prefix : prefixes){
				
				//counter++;
					
				extractor.reset(); //flushes motif map before new prefix is procesed
				//System.out.println("p: "+prefix);
			
				if (countNumberOfDegs(prefix)>motifSearchSpace.getMaxNumberOfDegeneratePositions()){
					continue prefixLoop;
				}
												
				for (int i=0; i<sequences.get(0).length(); i++){
					launchDiscoveryForPosition(i, factory, prefix, sequencesPreprocessed,forwardIndexes.get(i));
				}
				
				for (int i=0; i<revSequences.get(0).length(); i++){
					launchDiscoveryForPosition(i, factory, prefix, revSequencesPreprocessed,reverseIndexes.get(i));
				}
				
				
				/*if (counter%100==0){
					System.err.println("Prefix: " + counter + " "+(System.nanoTime()-start)/1000000);
					start = System.nanoTime();
				}*/
				
				

			}
		
		}
	}
	
	
	
	private ArrayList<IndexStructure> generateIndexes(ArrayList<Sequence> prepSeqs) {
		
		int numberOfInds=prepSeqs.get(0).length();
		ArrayList<IndexStructure> list = new ArrayList<IndexStructure>(numberOfInds);
		
		for (int i=0; i<numberOfInds; i++){
			
			List<Suffix> suffixes = new ArrayList<Suffix>();			
		
			for (int j=0; j<prepSeqs.size(); j++){
				suffixes.add(new Suffix(j,i));
			}
			
			list.add(iSFac.createIndexStructureForSuffixes(prepSeqs, suffixes));
		}
		return list;
		
	}

	private void launchDiscoveryForPosition(int i, MotifFactory factory, String prefix, ArrayList<Sequence> prepSeqs, IndexStructure ind){
				
		ISMonkey monkey = ind.getExactISMonkey(factory,motifSearchSpace.getMaximumDegeneracy());
		
		for (int p=0; p<prefix.length(); p++){
			monkey.jumpTo(prefix.charAt(p));
			if (!monkey.hasMatches()){
				return;
			}
		}
		
		//System.out.println("exploreSubtree for pos "+i+" prefix: ");
		exploreSubtreeFast(monkey);
	}

	private int countNumberOfDegs(String prefix) {
		
		Alphabet alph = motifSearchSpace.getAlphabet(); 
		int nDegs = 0;
		for (int i=0; i<prefix.length(); i++){
			nDegs+=alph.isDegenerate(prefix.charAt(i))?1:0;
		}
		return nDegs;
	}

	private SortedSet<String> generateAllPrefixes() {
		
		SortedSet<String> oldPrefs = new TreeSet<String>();
		SortedSet<String> newPrefs = new TreeSet<String>();
		SortedSet<String> dummySet;
		
		
		oldPrefs.add("");
		
		Alphabet alph = motifSearchSpace.getAlphabet();
		
		for (int i=0; i<prefixLength; i++){
			
			for (Character c : alph){
		
				for (String s : oldPrefs){
					newPrefs.add(s+c);
				}
				
			}
			
			dummySet= oldPrefs;
			oldPrefs = newPrefs;
			newPrefs = dummySet;
			newPrefs.clear();
				
				
		}
		
		return oldPrefs;
		
	}

	@Override
	public void setDataStructure(IndexStructureFactory indexStructureFactory) {
		
		this.iSFac = indexStructureFactory;
		
		
		
	}

	
}
