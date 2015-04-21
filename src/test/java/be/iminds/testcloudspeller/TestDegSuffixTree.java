package be.iminds.testcloudspeller;

import static org.junit.Assert.*;

import java.util.ArrayList;

import be.iminds.cloudspeller.indexing.GeneralizedSuffixTree;
import be.iminds.cloudspeller.indexing.ISMonkey;
import be.iminds.cloudspeller.indexing.SeqIDDecorationFactory;
import be.iminds.cloudspeller.input.BaseSequence;
import be.iminds.cloudspeller.input.Sequence;

import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

public class TestDegSuffixTree {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	
	}

	@Before
	public void setUp() throws Exception {
	}
	
	@Test
	/**
	 * Visual testing of suffix tree with dotformat
	 * dot -Tpdf testTree.txt > testTree.pdf
	 */
	public void testConstruction(){
		
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();
		
		//banana -> cagaga
		sequences.add(new BaseSequence("CAGAGA"));
		
		//build full
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(sequences,false,10,new SeqIDDecorationFactory());
		System.out.println(gst); 
		assertTrue(true); //TEST 30/5 OK
		
		
		//depth restriction=3
		gst = new GeneralizedSuffixTree(sequences,false,3,new SeqIDDecorationFactory());
		System.out.println(gst); 
		assertTrue(true); //TEST 30/5 OK
		
		//test reverse complement
		gst = new GeneralizedSuffixTree(sequences,true,3,new SeqIDDecorationFactory());
		System.out.println(gst); 
		assertTrue(true);//TEST 30/5 OK
		
		//multiple sequences + check sequence IDs
		sequences.clear();
		
		sequences.add(new BaseSequence("ACC")); 
		sequences.add(new BaseSequence("TTT"));
		sequences.add(new BaseSequence("ATA"));
		sequences.add(new BaseSequence("GGGC"));
		
		gst = new GeneralizedSuffixTree(sequences,false,2,new SeqIDDecorationFactory());
		System.out.println(gst); 
		assertTrue(true);//TEST 30/5 OK

	}
	
	@Test
	public void testPatternMatches(){
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();

		sequences.add(new BaseSequence("ACC")); 
		sequences.add(new BaseSequence("TTT"));
		sequences.add(new BaseSequence("ATA"));
		sequences.add(new BaseSequence("GGGC"));
		
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(sequences,false,10,new SeqIDDecorationFactory());

		//gst.getExactISMonkey(new IUPACFactory(IUPACType.FULL));
		
		//NON degenerate patterns
		assertEquals(gst.matchExactPattern(new IUPACMotif("A")).size(),3);
		assertEquals(gst.matchExactPattern(new IUPACMotif("C")).size(),3);
		assertEquals(gst.matchExactPattern(new IUPACMotif("G")).size(),3);
		assertEquals(gst.matchExactPattern(new IUPACMotif("T")).size(),4);
		assertEquals(gst.matchExactPattern(new IUPACMotif("AC")).size(),1);
		assertEquals(gst.matchExactPattern(new IUPACMotif("GG")).size(),2);
		
		//degenerate patterns
		assertEquals(gst.matchExactPattern(new IUPACMotif("N")).size(),13);
		
		assertEquals(gst.matchExactPattern(new IUPACMotif("M")).size(),6);
		assertEquals(gst.matchExactPattern(new IUPACMotif("R")).size(),6);
		assertEquals(gst.matchExactPattern(new IUPACMotif("W")).size(),7);
		assertEquals(gst.matchExactPattern(new IUPACMotif("S")).size(),6);
		assertEquals(gst.matchExactPattern(new IUPACMotif("Y")).size(),7);
		assertEquals(gst.matchExactPattern(new IUPACMotif("K")).size(),7);
		
		assertEquals(gst.matchExactPattern(new IUPACMotif("NNN")).size(),5);
	}
	
	@Test
	public void testMonkey(){
		
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();

		
		sequences.add(new BaseSequence("ACC")); 
		sequences.add(new BaseSequence("TTT"));
		sequences.add(new BaseSequence("ATA"));
		sequences.add(new BaseSequence("GGGC"));
		
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(sequences,false,10,new SeqIDDecorationFactory());

		int maxDeg=64;
		
		
		ISMonkey monkey = gst.getExactISMonkey(new IUPACFactory(IUPACType.FULL),maxDeg);
		
		//testing internal nodeInfos / motif trails / suffix grabbing
		
		assertTrue(monkey.grabInternalNodeInfo().toString().equals("[1111]"));
		monkey.jumpTo('A');
		assertTrue(monkey.getMotifTrail().toString().equals("A"));
		assertTrue(monkey.grabInternalNodeInfo().toString().equals("[1010]"));
		assertEquals(monkey.grabSuffixes().size(),3);
		monkey.jumpTo('C');
		assertTrue(monkey.getMotifTrail().toString().equals("AC"));
		assertTrue(monkey.grabInternalNodeInfo().toString().equals("[1000]"));
		assertEquals(monkey.grabSuffixes().size(),1);
		
		
		//ISMonkey monkeyClone = monkey.createClone();
		
		monkey.jumpTo('T');
		assertTrue(monkey.getMotifTrail().toString().equals("ACT"));
		assertTrue(monkey.grabInternalNodeInfo().toString().equals("[]"));
		assertEquals(monkey.grabSuffixes().size(),0);
		
		monkey.backtrack();
		
		monkey.jumpTo('C');
		assertTrue(monkey.getMotifTrail().toString().equals("ACC"));
		assertTrue(monkey.grabInternalNodeInfo().toString().equals("[1000]"));
		assertEquals(monkey.grabSuffixes().size(),1);
	
		//degenerate chars?
		
		ISMonkey newmonkey = gst.getExactISMonkey(new IUPACFactory(IUPACType.FULL),maxDeg);
		newmonkey.jumpTo('N');
		assertEquals(newmonkey.grabSuffixes().size(),13);
	}
	


	

	/*@Test OLD CODE
	public void testMonkeyUpdate(){
	
		//test1
		
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();
		
		sequences.add(new BaseSequence("ACC")); 
		sequences.add(new BaseSequence("TTT"));
		sequences.add(new BaseSequence("ATA"));
		sequences.add(new BaseSequence("GGGC"));
		
		int numChars = 0;
		for (int i=0; i<sequences.size(); i++){
			numChars+=sequences.get(i).length();
		}
		
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(sequences,false,10,new SeqIDDecorationFactory());

		int maxDeg=64;
		
		
		ISMonkey monkey = gst.getExactISMonkey(new IUPACFactory(IUPACType.FULL),maxDeg);
		
		ISMonkey monkeyA= monkey.createClone();
		monkeyA.jumpTo('A');
		ISMonkey monkeyC= monkey.createClone();
		monkeyC.jumpTo('C');
		ISMonkey monkeyG= monkey.createClone();
		monkeyG.jumpTo('G');
		ISMonkey monkeyT= monkey.createClone();
		monkeyT.jumpTo('T');
		
		List<ISMonkey> exactMonkeys = new LinkedList<ISMonkey>();
		exactMonkeys.add(monkeyA);
		exactMonkeys.add(monkeyC);
		exactMonkeys.add(monkeyG);
		exactMonkeys.add(monkeyT);
		
		
		//ISMonkey degMonkey = monkey.createDegenerateMonkey('N', exactMonkeys);
		
		/*assertTrue(degMonkey.getMotifTrail().toString().equals("N"));
		assertTrue(degMonkey.grabInternalNodeInfo().toString().equals("[1010, 1001, 0001, 0110]"));
		assertEquals(degMonkey.grabSuffixes().size(),numChars);*/
		
		
		//test2 sequences without a certain character
		/*
		sequences = new ArrayList<Sequence>();
		
		sequences.add(new BaseSequence("AAAAAAAAAAAAAAAA")); 
		sequences.add(new BaseSequence("CCCCCCCCCCCCCCCC"));
		sequences.add(new BaseSequence("GGGGGGGGGGGGGGGG"));
		sequences.add(new BaseSequence("AAAAAAAA"));
		
		numChars = 0;
		for (int i=0; i<sequences.size(); i++){
			numChars+=sequences.get(i).length();
		}
		
		gst = new GeneralizedSuffixTree(sequences,false,10,new SeqIDDecorationFactory());

		
		monkey = gst.getExactISMonkey(new IUPACFactory(IUPACType.FULL),maxDeg);
		
		monkeyA= monkey.createClone();
		monkeyA.jumpTo('A');
		monkeyC= monkey.createClone();
		monkeyC.jumpTo('C');
		monkeyG= monkey.createClone();
		monkeyG.jumpTo('G');
		monkeyT= monkey.createClone();
		monkeyT.jumpTo('T');
		
		exactMonkeys = new LinkedList<ISMonkey>();
		exactMonkeys.add(monkeyA);
		exactMonkeys.add(monkeyC);
		exactMonkeys.add(monkeyG);
		exactMonkeys.add(monkeyT);
		degMonkey = monkey.createDegenerateMonkey('N', exactMonkeys);
		
		assertTrue(degMonkey.getMotifTrail().toString().equals("N"));
		assertTrue(degMonkey.grabInternalNodeInfo().toString().equals("[1001, 0100, 0010]"));
		assertEquals(degMonkey.grabSuffixes().size(),numChars);
		
		
	}*/

	@Test
	public void testPreprocessor(){

		ArrayList<Sequence> sequences = new ArrayList<Sequence>();
		
		sequences.add(new BaseSequence("ACGTGGGT")); 
		sequences.add(new BaseSequence("TTTGGTGAA")); 
		sequences.add(new BaseSequence("AACCGGTTC")); 
		
		boolean withReverseComplement=false;
		
		GeneralizedSuffixTree index = new GeneralizedSuffixTree(sequences, withReverseComplement, 100, new SeqIDDecorationFactory());
		
		ArrayList<Sequence> prepSeqs = index.getPreprocessedSequences();
		
		for (int i=0; i<sequences.size(); i++){
			assertTrue(prepSeqs.get(i).toString().endsWith(GeneralizedSuffixTree.SENTINEL));
		}
		
		withReverseComplement=true;
		
		index = new GeneralizedSuffixTree(sequences, withReverseComplement, 100, new SeqIDDecorationFactory());
		prepSeqs = index.getPreprocessedSequences();
			
		assertEquals(prepSeqs.size(),6);
		assertTrue(prepSeqs.get(3).toString().equals("ACCCACGTs"));
		assertTrue(prepSeqs.get(4).toString().equals("TTCACCAAAs"));
		assertTrue(prepSeqs.get(5).toString().equals("GAACCGGTTs"));
		
		//check non alphabet char in sequences (has to replace with sentinels)
		sequences.add(new BaseSequence("ABCDEFG"));
		index = new GeneralizedSuffixTree(sequences, false, 100, new SeqIDDecorationFactory());
		prepSeqs = index.getPreprocessedSequences();
		assertTrue(prepSeqs.get(3).toString().equals("AsCsssGs"));
		
	}
	
	
		
	
	
	
	
}
