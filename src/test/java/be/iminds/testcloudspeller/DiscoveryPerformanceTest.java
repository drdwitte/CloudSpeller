package be.iminds.testcloudspeller;

import be.iminds.cloudspeller.indexing.BitSetDecorationFactory;
import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.NodeDecorationFactory;
import be.iminds.cloudspeller.input.BaseSequence;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;

import be.iminds.cloudspeller.motifalgorithms.DeNovoExactDiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.DiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.MotifContainer;
import be.iminds.cloudspeller.motifalgorithms.MotifSearchSpace;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;

import org.apache.commons.lang.NotImplementedException;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;

public class DiscoveryPerformanceTest {

	public static final String baseChars="ACGT";
	public static final boolean withReverseComplements = true;

	private static GeneFamily family;
	private static NodeDecorationFactory nodeDecoFac;
	
	
	public static void testConstruction(){
		throw new NotImplementedException();
	}
	
	public static void main(String [] args) {
		
		BLS.initializeBLSConstants(40,10,6);
		
		nodeDecoFac= new BitSetDecorationFactory();

		
		/*family = generateRandomDataset(500);
		nodeDecoFac= new BitSetDecorationFactory();
		System.out.println("DC ALPHABET k=12 seqLength=500, BitDecoration");
		checkLengthAndDegeneracyRange(new IUPACAlphabet(IUPACType.DONTCARES),12,12,64,64,4);
		System.out.println("");
		System.out.println("");
		*/
		
		family = generateRandomDataset(2000);
		System.out.println("DC ALPHABET k=12 seqLength=2000, BitDecoration");
		checkLengthAndDegeneracyRange(new IUPACAlphabet(IUPACType.DONTCARES),12,12,3,3);
		System.out.println("");
		System.out.println("");
		
		/*
		family = generateRandomDataset(500);
		System.out.println("TWOFOLDSANDN ALPHABET k=12 seqLength=500, BitDecoration");
		System.out.println("");
		checkLengthAndDegeneracyRange(new IUPACAlphabet(IUPACType.TWOFOLDSANDN),12,12,8,8,2);
		System.out.println("");
		System.out.println("");
		
		*/
		
		family = generateRandomDataset(2000);
		nodeDecoFac= new BitSetDecorationFactory();
		System.out.println("TWOFOLDSANDN ALPHABET k=6->12 seqLength=2000, BitDecoration");
		System.out.println("");
		checkLengthAndDegeneracyRange(new IUPACAlphabet(IUPACType.TWOFOLDSANDN),12,12,3,3);
		System.out.println("");
		System.out.println("");
		
	}
	
	
	
	
	public static void checkLengthAndDegeneracyRange(Alphabet alph, int kmin, int kmax, int minDegPos, int maxDegPos){
		for (int degPos=minDegPos; degPos<=maxDegPos; degPos++){
			System.out.println(((IUPACAlphabet)alph).getDegChars());
			for (int k=kmin; k<=kmax; k++){
				System.out.println("k= "+k+"\t #degPos="+degPos);
				double time = testRandomDataDiscovery(new MotifSearchSpace(kmin,k,degPos,alph));
				System.out.println(time+" milisec");
			}
		}
	}

	/**
	 * Test random data to estimate motif discovery on monocots (4 species, with rev comp)
	 * 
	 */
	public static double testRandomDataDiscovery(MotifSearchSpace searchSpace) {
		
		DiscoveryAlgorithm alg = new DeNovoExactDiscoveryAlgorithm(family.getSequences());
		BLSCalculator calculator=new BLSCalculator(family);
		calculator.setCutoff(new BLS(BLS.MIN));
		alg.setConservationScoreCalculator(calculator);
		int maxDepth=searchSpace.getMaxLength();
		alg.setDataStructure(new GSTFactory(maxDepth, withReverseComplements, nodeDecoFac));
		alg.setSearchSpace(searchSpace);
		long start=System.nanoTime();
		MotifContainer extractor = new MotifContainer();
		alg.setMotifExtractor(extractor);
		
		alg.runDiscovery(new IUPACFactory(IUPACType.FULL));
		long stop=System.nanoTime();
		System.out.println(((DeNovoExactDiscoveryAlgorithm)alg).getTempNumberOfMotifs()+ " motifs found");
		return (stop-start)/1000000;
	
	}

	private static GeneFamily generateRandomDataset(int sequenceLength) {
		GeneFamily family = new GeneFamily("TestFamily");
		family.setNewick("((g1:0.125,g2:0.125):0.25,(g3:0.125,g4:0.125):0.25);");
		
		family.addGeneSeq(new Gene("g1","kip"), new BaseSequence(generateRandomBaseString(sequenceLength)));
		family.addGeneSeq(new Gene("g2","konijn"), new BaseSequence(generateRandomBaseString(sequenceLength)));
		family.addGeneSeq(new Gene("g3","kameel"), new BaseSequence(generateRandomBaseString(sequenceLength)));
		family.addGeneSeq(new Gene("g4","kiwi"), new BaseSequence(generateRandomBaseString(sequenceLength)));

		return family;
	}

	private static String generateRandomBaseString(int sequenceLength) {
		StringBuffer s = new StringBuffer();
		
		for (int i=0; i<sequenceLength; i++){
			s.append(generateRandomBase());
		}
		
		return s.toString();
	}

	private static char generateRandomBase() {
		int index=(int)(Math.random()*4);
		return baseChars.charAt(index);
	}
}
