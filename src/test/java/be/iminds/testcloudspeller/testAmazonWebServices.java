package be.iminds.testcloudspeller;

import static org.junit.Assert.*;
import be.iminds.cloudspeller.indexing.BitSetDecoration;
import be.iminds.cloudspeller.indexing.BitSetDecorationFactory;
import be.iminds.cloudspeller.indexing.GeneralizedSuffixTree;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.indexing.NodeDecorationFactory;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.GeneFamily;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import be.iminds.cloudspeller.motifalgorithms.MotifSearchSpace;
import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;


import org.junit.BeforeClass;
import org.junit.Test;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.output.BLSConfGraphFactory;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;

/**
 * Test be.iminds.cloudspeller.output of CloudSpeller from AWS run on 4 gene families with k=5-6 on TWOFOLDSANDN
 * Here we try and reproduce the results using pattern matching and manual generation of
 * background models.
 * @author ddewitte
 *
 */

//extra test: aantal motieven uit elke gf:
// 197652 - 226335 - 255055 - 228903

public class testAmazonWebServices {

	private static ArrayList<GeneFamily> datasets = new ArrayList<GeneFamily>();
	private static Map<Motif,FreqVec> awsOutput = new HashMap<Motif,FreqVec>();	
	private static ArrayList<GeneralizedSuffixTree> indices = new ArrayList<GeneralizedSuffixTree>();
	private static ArrayList<BLSCalculator> blsCalculators =  new ArrayList<BLSCalculator>();
	private static MotifSearchSpace space = new MotifSearchSpace(5,6,2, new IUPACAlphabet(IUPACType.TWOFOLDSANDN));
	private static MotifFactory motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
	
	@BeforeClass
	public static void testSetup() {
		String dir = "AWSTest/";
		String [] files = {"iORTHO000100.txt", "iORTHO000101.txt", "iORTHO000102.txt", "iORTHO000104.txt"};
		BufferedReader in;
		System.out.print("Setup: generate families ...");
		for (String s : files){
			String filename = dir + s;
			try {
				in = new BufferedReader(new FileReader(filename));
				GeneFamily newFam = new GeneFamily(in); 
				datasets.add(newFam);
				
				in.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		System.out.println("done");
		
		//set parameters (see testcloudspeller)
		FreqVec.setNumberOfIntervals(6);
		BLS.initializeBLSConstants(40,10,6);
		boolean withReverseComplements=true;
		NodeDecorationFactory nodeDecorationFactory = new BitSetDecorationFactory();
		
		//build 4 gsts
		System.out.print("Building indices...");
		for (int i=0; i<datasets.size(); i++){
			indices.add(new GeneralizedSuffixTree(datasets.get(i).getSequences(), withReverseComplements, 15, nodeDecorationFactory));
			
		}
		System.out.println("done");
		
		//bls calculators
		System.out.print("Building bls machinery...");
		for (int i=0; i<datasets.size(); i++){
			blsCalculators.add(new BLSCalculator(datasets.get(i)));
		}
		System.out.println("done");
		
		//read aws be.iminds.cloudspeller.output
		System.out.print("Reading AWS be.iminds.cloudspeller.output...");
		String outputFilename = "outputTestSetAWS19juni.txt";
		try {
			in = new BufferedReader(new FileReader(outputFilename));
			BLSConfGraphFactory factory = new BLSConfGraphFactory();
			
			factory.getFreqVecsFromBuffer(in, awsOutput);
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("done");
		System.out.println("be.iminds.cloudspeller.output size: "+awsOutput.size());
		
		
	}
	
	@Test
	public void testOccurrences(){
		//generate N random patterns and perform double check
		int numberOfRandomPatterns = 10000;
		
		ArrayList<IUPACMotif> randomPatterns = generateRandomPatterns(numberOfRandomPatterns);
		
		for (IUPACMotif m : randomPatterns){
			FreqVec vAWS = testMotifAWS(m);
			FreqVec vPM = testMotifPM(m);
			/*System.out.println(m);
			System.out.println(vAWS);
			System.out.println(vPM);*/
			assertTrue(vAWS.equals(vPM));
		}
	}
	
	private FreqVec testMotifPM(IUPACMotif m) {
		
		
		FreqVec [] vecs = new FreqVec[datasets.size()]; 
		
		for (int i=0; i<indices.size(); i++){
			List<Suffix> suffixes =indices.get(i).matchExactPattern(m);
			
			if (suffixes==null){
				vecs[i]=new FreqVec();
				continue;
			}
			
			
			NodeDecoration b = new BitSetDecoration(datasets.get(i).getNumberOfGenes());
			b.processSuffixes(suffixes);
			BLSCalculator calculator = blsCalculators.get(i);
		
			BLS score = (BLS) calculator.calculateScore(b); 
			if (score!=null){
				vecs[i]=score.createFrequencyVector();
			} else {
				vecs[i]=new FreqVec();
			}
		}
		
		FreqVec total = new FreqVec();
		for (int i=0; i<vecs.length; i++){
			total.add(vecs[i]);
		}
		
		return total;
		
	}

	private FreqVec testMotifAWS(Motif m) {
		FreqVec res = awsOutput.get(m);
		if (res!=null){
			return res;
		} else {
			Motif complement = ((IUPACFactory)motifFactory).generateComplementMotif(m);
			res = awsOutput.get(complement);
			if (res!=null){
				return res;
			}
		}
		return new FreqVec();
	}

	private ArrayList<IUPACMotif> generateRandomPatterns(int numberOfRandomPatterns) {
		ArrayList<IUPACMotif> motifs = new ArrayList<IUPACMotif>();
		for (int i=0; i<numberOfRandomPatterns; i++){
			motifs.add(generateRandomPattern());
		}
		return motifs;
	}
	
	private IUPACMotif generateRandomPattern() {
				
		int k = generateRandomLength();
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

	private int generateRandomLength() {
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

	private char generateRandomChar(String s) {
		int index = generateRandomIndex(s.length());
		return s.charAt(index);
	}

	private int generateRandomIndex(int length) {
		return (int)(Math.random()*length);
	}

	
}
