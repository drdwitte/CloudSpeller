package be.iminds.testcloudspeller;

import static org.junit.Assert.*;

import org.apache.hadoop.io.BytesWritable;
import org.junit.Test;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.MotifFreqVec;

public class TestByteRepresentations {

	
	private static String createBinaryStringForByte(byte b){
		
		String trimmedBinaryRep = Integer.toBinaryString(b & 0xFF); //trimmed binary rep
		String strAtCorrectLength = String.format("%8s",trimmedBinaryRep); //reserve 8 chars
		return strAtCorrectLength.replace(' ', '0');

	}
	
	@Test
	/**
	 * Conversions: //http://mistupid.com/computers/binaryconv.htm
	 */
	public void testFreqVecs(){
		FreqVec.setNumberOfIntervals(3);
		
		//main case -> small numbers
		
		//[20,12,5] -> differences [8,7,5] -> revert [5,7,8] 
		//-> bytes [00000101,00000111,00001000]
		FreqVec vec = new FreqVec();
		vec.setFreq(0,20); vec.setFreq(1,12); vec.setFreq(2,5);
		byte [] bytesActual = vec.createBytesRepresentation();
		String [] bytesExpected = {"00000101","00000111","00001000"};
		for (int i=0; i<bytesActual.length; i++){
			String actual = createBinaryStringForByte(bytesActual[i]);
			assertEquals(bytesExpected[i], actual);
		}
		
		FreqVec vecDeserialized = FreqVec.createFreqVecFromBytes(0,bytesActual);
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			assertEquals(vec.getFreq(i),vecDeserialized.getFreq(i));
		}
//		
//		//special case 1 [20000,15000,3000] -> differences [5000,12000,3000]
//		// -> revert [3000,12000,5000]
//		
//		// 3000 = 0010111 0111000 -> 10111000, 00010111 
//		// 12000 = 1011101 1100000 -> 11100000, 01011101 
//		// 5000 = 0100111 0001000 -> 10001000, 00100111
//		//-> bytes [10111000, 00010111, 11100000, 01011101, 10001000, 00100111]
		
		vec = new FreqVec();
		vec.setFreq(0,20000); vec.setFreq(1,15000); vec.setFreq(2,3000);
		byte [] bytesActualSC1 = vec.createBytesRepresentation();;
		String [] bytesExpectedSC1 = {"10111000", "00010111", "11100000", "01011101", 
										"10001000", "00100111"};
		for (int i=0; i<bytesActualSC1.length; i++){
			String actual = createBinaryStringForByte(bytesActualSC1[i]);
			assertEquals(bytesExpectedSC1[i], actual);
		}
		
		vecDeserialized = FreqVec.createFreqVecFromBytes(0,bytesActualSC1);
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			assertEquals(vec.getFreq(i),vecDeserialized.getFreq(i));
		}
		
		//test shift
		int shift=2;
		byte [] bytesActualCopyShift = new byte[bytesActualSC1.length+shift];
		for (int i=0; i<bytesActualSC1.length; i++){
			bytesActualCopyShift[i+2]=bytesActualSC1[i];
		}
		
		for (int i=0; i<shift; i++){
			bytesActualCopyShift[i]=5;
		}
		
		vecDeserialized = FreqVec.createFreqVecFromBytes(shift,bytesActualCopyShift);
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			assertEquals(vec.getFreq(i),vecDeserialized.getFreq(i));
		}
//		
//		//special case 2 [18005,18000,2] -> differences [5,17798,2]
//		// -> revert [2,17798,5]
//		// 2 = 0000010 -> 00000010  
//		// 17998 = 0000001 0001100 1001110 -> 11001110, 10001100, 00000001 
		//  5 = 0000101-> 00000101 
//		//-> bytes [00000010, 11001110, 10001100, 00000001, 00000101]

		vec = new FreqVec();
		vec.setFreq(0,18005); vec.setFreq(1,18000); vec.setFreq(2,2);
		byte [] bytesActualSC2 = vec.createBytesRepresentation();
		String [] bytesExpectedSC2 = {"00000010", "11001110", "10001100", 
										"00000001", "00000101"};
		for (int i=0; i<bytesActualSC2.length; i++){
			String actual = createBinaryStringForByte(bytesActualSC2[i]);
			assertEquals(bytesExpectedSC2[i], actual);
		}
		
		vecDeserialized = FreqVec.createFreqVecFromBytes(0,bytesActualSC2);
		
		for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
			assertEquals(vec.getFreq(i),vecDeserialized.getFreq(i));
		}
		
		
		

		
	}
	
	@Test
	public void testMotifs(){
		
		//test compression for 2 different lengths
		

		//Test basepairs with different lengths
		
		IUPACFactory factory = new IUPACFactory(IUPACType.BASEPAIRS);
		factory.setMaxLength(12);
		IUPACMotif [] motifs12 =	{ new IUPACMotif("TATA"),
									new IUPACMotif("TTAACCGGT"),
									new IUPACMotif("AAAATTTTCCCC"),	
									new IUPACMotif("AAATTTTCCCC"),	
									(IUPACMotif)factory.createRandomMotif(5), 
									(IUPACMotif)factory.createRandomMotif(6), 
									(IUPACMotif)factory.createRandomMotif(7), 
									(IUPACMotif)factory.createRandomMotif(8), 
									(IUPACMotif)factory.createRandomMotif(9), 
									(IUPACMotif)factory.createRandomMotif(10), 
									(IUPACMotif)factory.createRandomMotif(11),
									(IUPACMotif)factory.createRandomMotif(12) 
									};
		
					
		for (int i=0; i<motifs12.length; i++){
			byte [] byteRep = factory.createBytesRepresentation(motifs12[i]);
			IUPACMotif m = (IUPACMotif) factory.createMotifFromBytes(0,byteRep);
			assertTrue(m.equals(motifs12[i]));
			assertEquals(byteRep.length,6);
		}
		
		//test degenerate motifs with different maxLength
		factory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
		factory.setMaxLength(15);
		IUPACMotif [] motifs15 = 	{ new IUPACMotif("TANAGW"), 
				(IUPACMotif)factory.createRandomMotif(5), 
				(IUPACMotif)factory.createRandomMotif(6), 
				(IUPACMotif)factory.createRandomMotif(7), 
				(IUPACMotif)factory.createRandomMotif(8), 
				(IUPACMotif)factory.createRandomMotif(9), 
				(IUPACMotif)factory.createRandomMotif(10), 
				(IUPACMotif)factory.createRandomMotif(11),
				(IUPACMotif)factory.createRandomMotif(12) 
									};
		
				
		for (int i=0; i<motifs15.length; i++){
			byte [] byteRep = factory.createBytesRepresentation(motifs15[i]);
			IUPACMotif m = (IUPACMotif) factory.createMotifFromBytes(0,byteRep);
			assertTrue(m.equals(motifs15[i]));
			assertEquals(byteRep.length,8);
		}
	}
	
	@Test 
	public void testMotifFreqVecs(){
		FreqVec.setNumberOfIntervals(5);
		IUPACFactory fac = new IUPACFactory(IUPACType.TWOFOLDSANDN);
		fac.setMaxLength(13);
		MotifFreqVec.setMotifFactory(fac);
		MotifFreqVec mFreq = new MotifFreqVec();
		mFreq.setMotif((IUPACMotif)fac.createRandomMotif(12));
		FreqVec vec = new FreqVec();
		int [] v = {100,10,2,1,0}; 
		vec.set(v);
		mFreq.setVec(vec);
		
		BytesWritable bW = mFreq.createBytesRepresentation();
		MotifFreqVec mFreqActual = MotifFreqVec.createMotifFreqVecFromBytes(bW);
		
		assertEquals(mFreq, mFreqActual);
		
		
		
		FreqVec.setNumberOfIntervals(6);
		fac.setMaxLength(7);
		MotifFreqVec.setMotifFactory(fac);
		mFreq = new MotifFreqVec();
		mFreq.setMotif((IUPACMotif)fac.createRandomMotif(7));
		vec = new FreqVec();
		int [] v2 = {100000,5,4,3,2,1}; 
		vec.set(v2);
		mFreq.setVec(vec);
		
		bW = mFreq.createBytesRepresentation();
		mFreqActual = MotifFreqVec.createMotifFreqVecFromBytes(bW);
		
		assertEquals(mFreq, mFreqActual);
		
		

		
	}
	
	
	
}
