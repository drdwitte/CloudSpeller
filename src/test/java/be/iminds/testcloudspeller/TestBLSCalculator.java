package be.iminds.testcloudspeller;

import static org.junit.Assert.*;
import be.iminds.cloudspeller.indexing.BitSetDecoration;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.indexing.SequenceIDSet;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.BaseSequence;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;

import java.util.ArrayList;
import java.util.List;
import org.junit.BeforeClass;
import org.junit.Test;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;
import be.iminds.cloudspeller.phylogenetics.NewickParser;

public class TestBLSCalculator {

	private static GeneFamily gf;
	private static String expectedSeqIDString;
	private static BLSCalculator calculator;
	
	@BeforeClass
	public static void testSetup() {
		
		gf = new GeneFamily("TestFam");
		Gene g1 = new Gene("g1","");
		Gene g2 = new Gene("g2","");
		Gene g3 = new Gene("g3","");
		Gene g4 = new Gene("g4","");
		Gene g5 = new Gene("g5","");
		Gene g6 = new Gene("g6","");
		
		BaseSequence seq = new BaseSequence("ACGT");
				
		gf.addGeneSeq(g1, seq);
		gf.addGeneSeq(g2, seq);
		gf.addGeneSeq(g3, seq);
		gf.addGeneSeq(g4, seq);
		gf.addGeneSeq(g5, seq);
		gf.addGeneSeq(g6, seq);
		
		String newick="((g1:0.125,(g2:0.025,g3:0.025):0.075):0.20,(g4:0.125,g5:0.125):0.20,g6:0.10);";
		gf.setNewick(newick);
		
		expectedSeqIDString="((0:0.125,(1:0.025,2:0.025):0.075):0.20,(3:0.125,4:0.125):0.20,5:0.10);";

	
	}
	
	
	@Test
	public void testGeneIDSubstitution() {
		String expectedSeqIDString="((0:0.125,(1:0.025,2:0.025):0.075):0.20,(3:0.125,4:0.125):0.20,5:0.10);";
		
		assertEquals(expectedSeqIDString, BLSCalculator.generateModifiedNewick(gf));
	}
	
	@Test
	public void testSubtrees(){
		String modifiedNewick= expectedSeqIDString.substring(1,expectedSeqIDString.length()-2);
		List<String> subtrees= NewickParser.getSubtrees(modifiedNewick);
		
		String [] subtreesByHand = {"(0:0.125,(1:0.025,2:0.025):0.075):0.20","(3:0.125,4:0.125):0.20","5:0.10"};
		int i=0;
		for (String s : subtrees){
			assertEquals(subtreesByHand[i++],s);
		}
	}

	

	@Test
	public void testBLSScoreCalculation(){
		calculator = new BLSCalculator(gf);
				
		//single branches:
		testSingleBranch(0,0.0);
		testSingleBranch(1,0.0);
		testSingleBranch(2,0.0);
		testSingleBranch(3,0.0);
		testSingleBranch(4,0.0);
		testSingleBranch(5,0.0);
		
		testTwoBranches(0,1,22.5);
		testTwoBranches(0,2,22.5);
		testTwoBranches(0,3,65.0);
		testTwoBranches(0,4,65.0);
		testTwoBranches(0,5,42.5);
		testTwoBranches(1,2,5.0);
		testTwoBranches(1,3,62.5);
		testTwoBranches(1,4,62.5);
		testTwoBranches(1,5,40.0);
		testTwoBranches(2,3,62.5);
		testTwoBranches(2,4,62.5);
		testTwoBranches(2,5,40.0);
		testTwoBranches(3,4,25.0);
		testTwoBranches(3,5,42.5);
		testTwoBranches(4,5,42.5);
		
		testFull();
		
	}
	
	private void testFull() {
		
		List<Suffix> suffixes = new ArrayList<Suffix>();
		int randSeqPos=500;
		for (int i=0; i<6; i++){
			suffixes.add(new Suffix(i,randSeqPos));
		}
			
		NodeDecoration nodeDeco1 = new SequenceIDSet(gf.getNumberOfGenes());
		nodeDeco1.processSuffixes(suffixes);

		NodeDecoration nodeDeco2 = new BitSetDecoration(gf.getNumberOfGenes());
		nodeDeco2.processSuffixes(suffixes);
		
		double expected = 100.0;
		double actualSeq = ((BLS)calculator.calculateScore(nodeDeco1)).doubleValue();
		double actualBits =((BLS)calculator.calculateScore(nodeDeco2)).doubleValue(); 
		assertEquals(expected, actualSeq,3);
		assertEquals(expected, actualBits,3);
		
		
		assertEquals(0,calculator.calculateScore(nodeDeco1).compareTo(new BLS(100)));
		assertEquals(0,calculator.calculateScore(nodeDeco2).compareTo(new BLS(100)));
	}
	
	private void testTwoBranches(int seqID1, int seqID2, double realBL) {
		int bl = (int)Math.round(realBL);
		List<Suffix> suffixes = new ArrayList<Suffix>();
		int randSeqPos=500;
		suffixes.add(new Suffix(seqID1,randSeqPos));
		suffixes.add(new Suffix(seqID2,randSeqPos));
		
		NodeDecoration nodeDeco1 = new SequenceIDSet(gf.getNumberOfGenes());
		nodeDeco1.processSuffixes(suffixes);

		NodeDecoration nodeDeco2 = new BitSetDecoration(gf.getNumberOfGenes());
		nodeDeco2.processSuffixes(suffixes);
		
		
		double expected = bl;
		double actualSeq = ((BLS)calculator.calculateScore(nodeDeco1)).doubleValue();
		double actualBits =((BLS)calculator.calculateScore(nodeDeco2)).doubleValue(); 
	
		assertEquals(expected, actualSeq,1.0); //rounding errors may cause max diff of 1.0
		assertEquals(expected, actualBits,1.0);
		
		
		
	}

	public void testSingleBranch(int seqID,double realBL){
		int bl = (int)Math.round(realBL);
		List<Suffix> suffixes = new ArrayList<Suffix>();
		int randSeqPos=500;
		suffixes.add(new Suffix(seqID,randSeqPos));
		
		NodeDecoration nodeDeco1 = new SequenceIDSet(gf.getNumberOfGenes());
		nodeDeco1.processSuffixes(suffixes);
		
		NodeDecoration nodeDeco2 = new BitSetDecoration(gf.getNumberOfGenes());
		nodeDeco2.processSuffixes(suffixes);
		

		assertEquals(0,calculator.calculateScore(nodeDeco1).compareTo(new BLS(bl)));
		assertEquals(0,calculator.calculateScore(nodeDeco2).compareTo(new BLS(bl)));
	}
	
	

}
