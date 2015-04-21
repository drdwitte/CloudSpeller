package be.iminds.testcloudspeller;

import static org.junit.Assert.*;
import be.iminds.cloudspeller.indexing.BitSetDecoration;
import be.iminds.cloudspeller.indexing.BitSetDecorationFactory;
import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.GeneralizedSuffixTree;
import be.iminds.cloudspeller.indexing.ISMonkey;
import be.iminds.cloudspeller.indexing.IndexStructure;
import be.iminds.cloudspeller.indexing.IndexStructureFactory;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.indexing.NodeDecorationFactory;
import be.iminds.cloudspeller.indexing.SeqIDDecorationFactory;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.BaseSequence;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import be.iminds.cloudspeller.motifalgorithms.ABMotifContainer;
import be.iminds.cloudspeller.motifalgorithms.AlBasedExactDiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.DeNovoExactDiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.MotifContainer;
import be.iminds.cloudspeller.motifalgorithms.MotifSearchSpace;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import org.junit.Test;

import be.iminds.cloudspeller.driver.DistributedABPatternMatcher;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;
import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculator;

public class TestExactDiscovery {
	
	@Test
	public void testExactDiscovery(){
		
		
		GeneFamily family = new GeneFamily("TestFamily");
		family.addGeneSeq(new Gene("g1","zma"), new BaseSequence("AACGTTCA"));
		family.addGeneSeq(new Gene("g2","osa"), new BaseSequence("TTAAGCCC"));
		family.addGeneSeq(new Gene("g3","bdi"), new BaseSequence("TTCAAGCT"));
		family.addGeneSeq(new Gene("g4","sbi"), new BaseSequence("TTAAGCCTT"));
		
		family.setNewick("((g1:0.086,g4:0.086):0.2366,(g2:0.2688,g3:0.2688):0.0538);"); //newick from monocots paper
		
		DeNovoExactDiscoveryAlgorithm discoveryAlgorithm = new DeNovoExactDiscoveryAlgorithm(family.getSequences());
		
		ConservationScoreCalculator c = new BLSCalculator(family);
		c.setCutoff(new BLS(70));
		discoveryAlgorithm.setConservationScoreCalculator(c);
		
		NodeDecorationFactory fac = new SeqIDDecorationFactory();
		IndexStructureFactory indexStructureFactory = new GSTFactory(10,false,fac);
		discoveryAlgorithm.setDataStructure(indexStructureFactory);
		
		Alphabet motifAlphabet = new IUPACAlphabet(IUPACType.BASEPAIRS);
		int kmin=4; int kmax=4; int degPos=0;
		MotifSearchSpace motifSearchSpace = new MotifSearchSpace(kmin,kmax,degPos,motifAlphabet);
		discoveryAlgorithm.setSearchSpace(motifSearchSpace);
		
		MotifContainer extractor = new MotifContainer();
		discoveryAlgorithm.setMotifExtractor(extractor);
		discoveryAlgorithm.runDiscovery(new IUPACFactory(IUPACType.BASEPAIRS));
		Map<Motif,ConservationScore> motifMap = extractor.getMotifMap();
		
		
		//all 4mers (by hand)
		
		//AACG: {0}
		//AAGC: {1,2,3} -> 100 - 8.6 = 91.4 -> 91%
		//ACGT: {0}
		//AGCC: {1,3}	-> 100 - 26.88 - 8.6 = 100 - 35.48 = 64.52 -> 65%
		//AGCT: {2} 
		//CAAG: {2} 
		//CCTT: {3}
		//CGTT: {0}
		//GCCC: {1} 
		//GCCT: {3}
		//GTTC: {0}
		//TAAG: {1,3}	-> 65%
		//TTAA: {1,3}	-> 65%
		//TTCA: {0,2}	-> 65%
		//TCAA: {2} 
		
		assertTrue(motifMap.get(new IUPACMotif("AAGC")).equals(new BLS(91)));
	
		c.setCutoff(new BLS(60));
		
		discoveryAlgorithm.runDiscovery(new IUPACFactory(IUPACType.BASEPAIRS));
		motifMap = extractor.getMotifMap();
		
		
		assertTrue(motifMap.get(new IUPACMotif("AGCC")).equals(new BLS(65)));
		assertTrue(motifMap.get(new IUPACMotif("TAAG")).equals(new BLS(65)));
		assertTrue(motifMap.get(new IUPACMotif("TTAA")).equals(new BLS(65)));
		assertTrue(motifMap.get(new IUPACMotif("TTCA")).equals(new BLS(65)));
		
	
		motifAlphabet = new IUPACAlphabet(IUPACType.DONTCARES);
		degPos=1;
		motifSearchSpace = new MotifSearchSpace(kmin,kmax,degPos,motifAlphabet);
		discoveryAlgorithm.setSearchSpace(motifSearchSpace);
		c.setCutoff(new BLS(100));
		
		discoveryAlgorithm.runDiscovery(new IUPACFactory(IUPACType.DONTCARES));
		motifMap = extractor.getMotifMap();
		
		IndexStructure index =indexStructureFactory.createIndexStructure(family.getSequences());
		
	    
	    assertEquals(index.matchExactPattern(new IUPACMotif("TNAA")).size(),3);
	    assertEquals(index.matchExactPattern(new IUPACMotif("NCAA")).size(),1);
	    assertEquals(index.matchExactPattern(new IUPACMotif("NTAA")).size(),2);
	    assertEquals(index.matchExactPattern(new IUPACMotif("TTNA")).size(),4);
	    
	    int maxDeg=64;
	    
	    ISMonkey monkey = index.getExactISMonkey(new IUPACFactory(IUPACType.DONTCARES),maxDeg);
	    monkey.jumpTo('T');
	    assertTrue(monkey.grabInternalNodeInfo().toString().equals("[1111]"));
	    
	    
	    monkey.jumpTo('N');
	    assertTrue(monkey.grabInternalNodeInfo().toString().equals("[0101, 1010, 1111]"));
	    
	    monkey.jumpTo('A');
	    assertTrue(monkey.grabInternalNodeInfo().toString().equals("[0101, 1010, 0101]"));
	    
	    monkey.jumpTo('A');
	    assertTrue(monkey.grabInternalNodeInfo().toString().equals("[0010, 0101]"));
	    
    
	}
	
	//NOTE with deg=3 pos / k=6-7 takes approx 20 sec
	@Test
	public void testCompareDiscoveryWithPatternMatchingOnRandomData(){
		
		MotifFactory motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
		
		int numberOfSequences=4;
		int seqLength=500;
		int maxTreeDepth=10;
		boolean withReverseComplements=true;
			
		//de novo discovery approach 
		
		GeneFamily family = new GeneFamily("TestFamily");
		family.addGeneSeq(new Gene("g1","zma"), new BaseSequence(generateRandomDNASequence(seqLength)));
		family.addGeneSeq(new Gene("g2","osa"), new BaseSequence(generateRandomDNASequence(seqLength)));
		family.addGeneSeq(new Gene("g3","bdi"), new BaseSequence(generateRandomDNASequence(seqLength)));
		family.addGeneSeq(new Gene("g4","sbi"), new BaseSequence(generateRandomDNASequence(seqLength)));
		
		family.setNewick("((g1:0.086,g4:0.086):0.2366,(g2:0.2688,g3:0.2688):0.0538);"); //newick from monocots paper
		
		//System.out.println(family);
		
		
		DeNovoExactDiscoveryAlgorithm discoveryAlgorithm = new DeNovoExactDiscoveryAlgorithm(family.getSequences());
		
		BLSCalculator calculator = new BLSCalculator(family);
		int blsCutoff=40;
		calculator.setCutoff(new BLS(blsCutoff));
		discoveryAlgorithm.setConservationScoreCalculator(calculator);
		
		NodeDecorationFactory fac = new BitSetDecorationFactory();//new BitSetDecorationFactory();
		IndexStructureFactory indexStructureFactory = new GSTFactory(maxTreeDepth,withReverseComplements,fac);
		discoveryAlgorithm.setDataStructure(indexStructureFactory);
		
		Alphabet motifAlphabet = new IUPACAlphabet(IUPACType.TWOFOLDSANDN);
		int kmin=6; int kmax=7; int degPos=2;
		MotifSearchSpace motifSearchSpace = new MotifSearchSpace(kmin,kmax,degPos,motifAlphabet);
		discoveryAlgorithm.setSearchSpace(motifSearchSpace);
		
		MotifContainer extractor = new MotifContainer();
		discoveryAlgorithm.setMotifExtractor(extractor);
		discoveryAlgorithm.runDiscovery(new IUPACFactory(IUPACType.TWOFOLDSANDN));
		Map<Motif,ConservationScore> motifMap = extractor.getMotifMap();
		
		//pattern matching approach
		
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(family.getSequences(),withReverseComplements,maxTreeDepth,new SeqIDDecorationFactory());
		
		NodeDecoration deco = new SeqIDDecorationFactory().createNodeDecoration(numberOfSequences);//fac.createNodeDecoration(numberOfSequences);
		
		for (Map.Entry<Motif,ConservationScore> e : motifMap.entrySet()){
			IUPACMotif pattern = (IUPACMotif) e.getKey();
			ConservationScore deNovoScore = e.getValue();
			List<Suffix> suffixes = gst.matchExactPattern(pattern);
			deco.processSuffixes(suffixes);
			ConservationScore pmScore = calculator.calculateScore(deco);
			
			
			if (!pmScore.equals(deNovoScore)){
				System.out.println(pattern);
				System.out.println(pmScore+"\t"+deNovoScore);
				System.out.println(suffixes);
				System.out.println(deco);
			}
			assertTrue(pmScore.equals(deNovoScore));
		}
		
		
		//also check random patterns:
		
		List<IUPACMotif> patterns = generateRandomPatterns(1000,motifSearchSpace, motifFactory);
		for (IUPACMotif p : patterns){
			List<Suffix> suffixes = gst.matchExactPattern(p);
			ConservationScore discoveryScore = motifMap.get((Motif)p);
						
			if (suffixes==null){ //score = 0
				assertNull(discoveryScore);
			
			} else { //calculate score
				NodeDecoration b = new BitSetDecoration(family.getNumberOfGenes());
				b.processSuffixes(suffixes);
				BLS score = (BLS) calculator.calculateScore(b); 

				if (score!=null){
					assertNotNull(discoveryScore);
					assertEquals(0,score.compareTo(discoveryScore));
				} else {
					assertNull(discoveryScore);
				}
			}
		}

	}
	
	private ArrayList<IUPACMotif> generateRandomPatterns(int numberOfRandomPatterns, MotifSearchSpace space, MotifFactory motifFactory) {
		ArrayList<IUPACMotif> motifs = new ArrayList<IUPACMotif>();
		for (int i=0; i<numberOfRandomPatterns; i++){
			motifs.add(generateRandomPattern(space,motifFactory));
		}
		return motifs;
	}
	
	private IUPACMotif generateRandomPattern(MotifSearchSpace space, MotifFactory motifFactory) {
		
		int k = generateRandomLength(space);
		int degPos = space.getMaxNumberOfDegeneratePositions();
		int kExact = k - degPos;
		
		String exactChars = generateRandomString("ACGT",kExact);
		StringBuilder sb = new StringBuilder(exactChars);
		
		for (int i=0; i<degPos; i++){
			char c = generateRandomChar(space.getAlphabet().getAllChars());
			int randPos = generateRandomIndex(sb.length());
			sb.insert(randPos, c);
		}
		
		return (IUPACMotif)motifFactory.createMotifFromString(sb.toString());
		
	}
	
	private int generateRandomLength(MotifSearchSpace space) {
		int min = space.getMinLength();
		int numPos= space.getMaxLength()-min+1;
		int offset = (int)(Math.random()*numPos);
		return min+offset;
	}

	private String generateRandomString(String alphabet,
			int length) {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<length; i++){
			sb.append(generateRandomChar(alphabet));
		}
		return sb.toString();
	}

	private int generateRandomIndex(int length) {
		return (int)(Math.random()*length);
	}

	
	public static String generateRandomDNASequence(int length){
		StringBuffer buffer = new StringBuffer();
		for (int i=0; i<length; i++){
			buffer.append(generateRandomChar("ACGT"));
		}
		return buffer.toString();
	}
	
	public static char generateRandomChar(String alphabet){
		  int index=(int)Math.floor(Math.random()*alphabet.length());
		  return alphabet.charAt(index);
	}
	
	
	@Test
	public void AlignmentBasedDiscoveryAndPM(){
		
		//MotifFactory motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
		
		int numberOfSequences=4;
		//int seqLength=500;
		int maxTreeDepth=10;
		boolean withReverseComplements=true;
		
		GeneFamily family = new GeneFamily("TestFamily");
		family.addGeneSeq(new Gene("g1","zma"), new BaseSequence("AC-TTATAACGTA"));
		family.addGeneSeq(new Gene("g2","osa"), new BaseSequence("CG-ATATATTATA"));
		family.addGeneSeq(new Gene("g3","bdi"), new BaseSequence("G-ACTATATTATA"));
		family.addGeneSeq(new Gene("g4","sbi"), new BaseSequence("-ACGTATATA-GT"));
		
		family.setNewick("((g1:0.086,g4:0.086):0.2366,(g2:0.2688,g3:0.2688):0.0538);"); //newick from monocots paper
		
		//System.out.println(family);
		
		System.out.println("BEGIN");
		
		DeNovoExactDiscoveryAlgorithm discoveryAlgorithm = new AlBasedExactDiscoveryAlgorithm(family.getSequences());
		
		BLSCalculator calculator = new BLSCalculator(family);
		int blsCutoff=40;
		calculator.setCutoff(new BLS(blsCutoff));
		discoveryAlgorithm.setConservationScoreCalculator(calculator);
		
		NodeDecorationFactory fac = new BitSetDecorationFactory();
		IndexStructureFactory indexStructureFactory = new GSTFactory(maxTreeDepth,withReverseComplements,fac);
		discoveryAlgorithm.setDataStructure(indexStructureFactory);
		
		Alphabet motifAlphabet = new IUPACAlphabet(IUPACType.TWOFOLDSANDN);
		int kmin=4; int kmax=6; int degPos=1;
		MotifSearchSpace motifSearchSpace = new MotifSearchSpace(kmin,kmax,degPos,motifAlphabet);
		discoveryAlgorithm.setSearchSpace(motifSearchSpace);
		
		ABMotifContainer extractor = new ABMotifContainer();
		discoveryAlgorithm.setMotifExtractor(extractor);
		discoveryAlgorithm.runDiscovery(new IUPACFactory(IUPACType.TWOFOLDSANDN));
		Map<Motif,ConservationScore> motifMap = extractor.getMotifMap();
		
		assertTrue(motifMap.get(new IUPACMotif("TATA")).equals(new BLS(100)));
		
		assertTrue(motifMap.get(new IUPACMotif("TATA")).equals(new BLS(100))); //{1234}
		assertTrue(motifMap.get(new IUPACMotif("TATAT")).equals(new BLS(91))); //{234}
		assertTrue(motifMap.get(new IUPACMotif("TATATT")).equals(new BLS(54))); //{23}
		
		//reverse comps
		assertTrue(motifMap.get(new IUPACMotif("TATA")).equals(new BLS(100))); // {}
		assertTrue(motifMap.get(new IUPACMotif("ATATA")).equals(new BLS(91))); // {}
		assertTrue(motifMap.get(new IUPACMotif("AATATA")).equals(new BLS(54))); // {}
		
		//deg variant
		assertTrue(motifMap.get(new IUPACMotif("RTATA")).equals(new BLS(91))); // {}

		
		//check if misaligned is not counted:
		

		assertNull(motifMap.get(new IUPACMotif("CGTA"))); // {} (occurs alfree but not aligned!

		
		
		
		//pattern matching approach
		
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(family.getSequences(),withReverseComplements,maxTreeDepth,new SeqIDDecorationFactory());
		
		NodeDecoration deco = new SeqIDDecorationFactory().createNodeDecoration(numberOfSequences);//fac.createNodeDecoration(numberOfSequences);
		
		for (Map.Entry<Motif,ConservationScore> e : motifMap.entrySet()){
			IUPACMotif pattern = (IUPACMotif) e.getKey();
			ConservationScore deNovoScore = e.getValue();
			
			List<Suffix> suffixes = gst.matchExactPattern(pattern);
			
			assertNotNull(suffixes);			
			
			Map<Integer,List<Suffix>> alignedSuff = DistributedABPatternMatcher.GFABPatternMatcher.getAlignedSuffixesBySequencePosition(suffixes);
			
			ConservationScore pmScore = new BLS(0);
			
			List<Suffix> alSuffixes = null;
			String decoBestStr = null;
			for (Map.Entry<Integer,List<Suffix>> entry : alignedSuff.entrySet()){
				/*System.out.println(entry.getKey());
				System.out.println(entry.getValue());*/
				deco.processSuffixes(entry.getValue());
				BLS score = (BLS) calculator.calculateScore(deco);
				if (score==null)
					continue;
				if (score.compareTo(pmScore)>0){ //in principle this is done with the ABMotifContainer
					pmScore = score;
					alSuffixes=entry.getValue();
					decoBestStr = deco.toString();
				}
				
			}
		
			if (!pmScore.equals(deNovoScore)){
				System.out.println(pattern);
				System.out.println(pmScore+"\t"+deNovoScore);
				System.out.println(alSuffixes);
				System.out.println(decoBestStr);
			}
			assertTrue(pmScore.equals(deNovoScore));
		}
		
		System.out.println("END");
		
	}
	
					
	
	
	
	
	
	
	

	
	
	
	
	
}
