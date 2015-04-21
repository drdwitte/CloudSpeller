package be.iminds.testcloudspeller;

import static org.junit.Assert.*;
import org.junit.Test;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.postprocessing.JacardDistanceCalculator;
import be.iminds.cloudspeller.postprocessing.MotifAligner;

public class testMotifAligner {

	@Test
	public void testMotifDistances(){
		
		Alphabet alphabet = new IUPACAlphabet(IUPACType.TWOFOLDSANDN);
		MotifAligner aligner = new MotifAligner(alphabet, new JacardDistanceCalculator(alphabet));
		
		
		double actual=0.0; double expected=0.0;
		
		//substring test distance = 3 x d(C,.) = 3 x d(C,N) = 3 x (1.0 - 1/4) = 9/4
		actual = aligner.calculateSimilarityDistance("ATA","ATACCC");
		expected = 9.0/4;
		assertEquals(expected, actual, 3);
		
		//test1 ingnore boundary Ns
		actual = aligner.calculateSimilarityDistance("ATA","ATANN");
		expected = 0.0;
		assertEquals(expected, actual, 3);
		
		//test2 ingnore boundary Ns
		actual = aligner.calculateSimilarityDistance("CTT","NCTTNN");
		expected = 0.0;
		assertEquals(expected, actual, 3);
		
		//test weak boundaries: 
		actual = aligner.calculateSimilarityDistance("GGC","WWGGCNR");
		expected = 3.0/2.0;
		assertEquals(expected, actual, 3);
		
		//test regular hamming distance: 1 mismatch
		actual = aligner.calculateSimilarityDistance("CATACGC","CATATGC");
		expected = 1.0;
		assertEquals(expected, actual, 3);
		
		//test regular hamming distance: degenerate + 1 mismatch
		actual = aligner.calculateSimilarityDistance("CNTACGC","CATATGC");
		expected = 1.75;
		assertEquals(expected, actual, 3);
		
		//one complex one for fun
		actual = aligner.calculateSimilarityDistance("CNTACGC","CATATGC");
		expected = 1.75;
		assertEquals(expected, actual, 3);
		
	}
}
