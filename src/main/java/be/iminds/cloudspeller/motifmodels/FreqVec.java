package be.iminds.cloudspeller.motifmodels;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.ListIterator;


public class FreqVec {
	
	private static final int maxFreqInByte = (1 << 7) -1;
	private static final byte continuationBit = (byte)(1 << 7);
	private static final String delimiter=",";
	private static int numInterv;

	private int [] vec;
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof FreqVec){
			
			FreqVec other = (FreqVec) obj;
			int [] valuesOther = other.vec;
			for (int i=0; i<vec.length; i++){
				if (vec[i]!=valuesOther[i]){
					return false;
				}
			}
			return true;
		} 
		return false;
		
	}

	
	/**
	 * Default constructor
	 */
	public FreqVec(){
		vec = new int[numInterv];
	}
	

	//SETTERS
	
	public static void setNumberOfIntervals(int nInterv){
		numInterv=nInterv;
	}
	
	public void setFreq(int i, int f){
		vec[i]=f;
	}
	
	public void set(int [] v){
		for (int i=0; i<v.length; i++){
			setFreq(i,v[i]);
		}
	}
	
	
	//GETTERS
	
	public int getFreq(int i){
		return vec[i];
	}
	
	public static int getNumberOfIntervals() {
		return numInterv;
	}
	
	//IO
	
	public static FreqVec createFreqVecFromString(String s){
		String [] stringFreqs=s.split(delimiter);
		FreqVec freqVec= new FreqVec();
		for (int i=0; i<stringFreqs.length; i++){
			int f=Integer.parseInt(stringFreqs[i]);
			freqVec.setFreq(i, f);
		}
		return freqVec;
	}
	
	public String createStringRepresentation(){
		StringBuilder s = new StringBuilder();
		for (int i=0; i<vec.length; i++){
			s.append(vec[i]);
			s.append(delimiter);
		}
		return s.toString();
	}
	
	public static FreqVec createFreqVecFromBytes(int offset, byte [] byteArray){
		
		FreqVec vec = new FreqVec();
		int bAIndex=offset;
		int previousFreq=0;
		for (int i=0; i<numInterv; i++){
					
			int freq = 0;
			int sevenBits = 0;
			int shiftDepth= 0;
			
			while(isIncomplete(byteArray[bAIndex])){
				sevenBits = byteArray[bAIndex++] & maxFreqInByte;
				freq = freq | (sevenBits << (shiftDepth++)*7);
			}
			
			sevenBits = byteArray[bAIndex++] & maxFreqInByte;
			freq = (freq | (sevenBits << shiftDepth*7)) + previousFreq;
			vec.setFreq(numInterv-i-1,freq);
			previousFreq=freq;
			
			
		}
		return vec;
	}
	
	
	private static boolean isIncomplete(byte b){
		return (b & 0x80)!=0;
	}
	
	/**
	 * Use variable length integer encoding -> bytes representation using continuation
	 * bit => only 7 bits to represent number:
	 * For the number we choose differences, this increases the probability that every
	 * number can be represented using 1 byte
	 * @return byte array representation of a frequency vector
	 * 
	 */

	public byte [] createBytesRepresentation(){
	
		LinkedList<Byte> l = new LinkedList<Byte>();
		addBytesRep(vec[vec.length-1],l);
		for (int i=vec.length-2; i>=0; i--){
			addBytesRep(vec[i]-vec[i+1],l);
			
		}
		
		byte [] bytesRep = new byte[l.size()];
		
		ListIterator<Byte> iterator = l.listIterator(); 
		int c=0;
		while (iterator.hasNext()){
			bytesRep[c++]=iterator.next();
		}
		return bytesRep;
	}
	
	private void addBytesRep(int integerValue, LinkedList<Byte> l) {

		int sevenBits = integerValue & maxFreqInByte; // value & 0111 1111
		int rem = (integerValue >> 7);
		while (rem!=0){
			l.add((byte)(continuationBit|sevenBits));
			sevenBits = rem & maxFreqInByte;
			rem >>=7;
		}
		l.add((byte)sevenBits); //add 0_sevenBits 
		return;
	}

	//METHODS
	
	public void add(FreqVec v){
		for (int i=0; i<numInterv; i++){
			vec[i]+=v.vec[i];
		}
	}


	@Override
	public String toString() {
		return Arrays.toString(vec);
	}


	//DEBUG
	/*private static String createBinaryStringForByte(byte b){
		
		String trimmedBinaryRep = Integer.toBinaryString(b & 0xFF); //trimmed binary rep
		String strAtCorrectLength = String.format("%8s",trimmedBinaryRep); //reserve 8 chars
		return strAtCorrectLength.replace(' ', '0');

	}*/

}
