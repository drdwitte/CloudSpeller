package be.iminds.cloudspeller.motifmodels;

import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

public class IUPACFactory extends MotifFactory {

	private Alphabet alphabet;
	private int numberOfBytesForMotif=6;
	private int maxLength=12;
	private static final char noChar = ' ';
	private CharPair [] byteDiMerMap;
	private Map<CharPair,Byte> dimerByteMap = null;
	
	
	public IUPACFactory (IUPACType type){
		alphabet = new IUPACAlphabet(type);
		initializeByteCompressionUtilities();
	}
	
	private void initializeByteCompressionUtilities() {
		
		dimerByteMap = new HashMap<CharPair,Byte>();
		int alphSize = alphabet.getAllChars().length();
		int mapSize = alphSize*alphSize + alphSize + 1;
		byteDiMerMap = new CharPair[mapSize];
		int value=0;	
				
		for (Character c1 : alphabet){
			for (Character c2 : alphabet){
				CharPair p = new CharPair(c1, c2);
				byteDiMerMap[value]=p;
				dimerByteMap.put(p,(byte)value++);
				
			}
			
			CharPair p = new CharPair(c1,noChar);
			byteDiMerMap[value]=p;
			dimerByteMap.put(p,(byte)value++);
		}
		
		
		CharPair p = new CharPair(noChar,noChar);
		byteDiMerMap[value]=p;
		dimerByteMap.put(new CharPair(noChar,noChar),(byte)value);
	}

	@Override
	public Motif createEmptyMotif() {
		return new IUPACMotif("");
	}

	@Override
	public Motif createMotifFromBytes(int offset, byte [] b) {
		
		StringBuilder sb = new StringBuilder();
		for (int i=offset; i<numberOfBytesForMotif; i++){
			sb.append(byteDiMerMap[byteToIndex(b[i])]);
		}
		
		return new IUPACMotif(sb.toString().trim());
	}
	
	private int byteToIndex(byte b){
		return b & 0xFF;
	}


	@Override
	public Motif createMotifFromString(String s) {
		return new IUPACMotif(s);
	}

	@Override
	public Alphabet getAlphabet() {
		return alphabet;
	}

	@Override
	public int getNumberOfBytesForMotif() {
		return numberOfBytesForMotif;
	}

	@Override
	public byte[] createBytesRepresentation(Motif m) {
		
		byte [] byteRep = new byte[numberOfBytesForMotif];
		int byteIndex=0;
		int i=0;
		for (; i<m.length()-1; i+=2){
			CharPair pair = new CharPair(m.charAt(i),m.charAt(i+1));
			byteRep[byteIndex++]=dimerByteMap.get(pair);
		}
		
		
		if (i==m.length()-1){
			CharPair pair = new CharPair(m.charAt(i),noChar);
			byteRep[byteIndex++]=dimerByteMap.get(pair);
		}
		
		//add empty chars
		byte noCharRep = dimerByteMap.get(new CharPair(noChar,noChar));
		for (; byteIndex<numberOfBytesForMotif; byteIndex++){
			byteRep[byteIndex]=noCharRep;
		}
		
		return byteRep;
	}

	@Override
	public String createStringRepresentation(Motif m) {
		return m.toString();
	}

	@Override
	public int getMaxLength() {
		return maxLength;
	}

	@Override
	public void setMaxLength(int max) {
		maxLength=max;
		numberOfBytesForMotif=max/2+max%2;
	}
	
	
	public Motif generateComplementMotif(Motif m){
		String mStr = m.toString();
		StringBuilder compBuilder = new StringBuilder();
		for (int i=mStr.length()-1; i>=0; i--){
			compBuilder.append(alphabet.getComplement(mStr.charAt(i)));
		}
		return new IUPACMotif(compBuilder.toString(),((IUPACMotif)m).numberOfDegPositions());
		
	}
	
	//DEBUG
	/*private static String createBinaryStringForByte(byte b){
		
		String trimmedBinaryRep = Integer.toBinaryString(b & 0xFF); //trimmed binary rep
		String strAtCorrectLength = String.format("%8s",trimmedBinaryRep); //reserve 8 chars
		return strAtCorrectLength.replace(' ', '0');

	}*/
}

class CharPair {
	
	private char c1;
	private char c2;
	
	public CharPair(char c1, char c2){
		this.c1=c1; this.c2=c2;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof CharPair){
			CharPair cp = (CharPair) obj;
			return this.c1 == cp.c1 && this.c2==cp.c2;
		}
		return false;
	}
	@Override
	public int hashCode() {
		return toString().hashCode();
	}
	@Override
	public String toString() {
		return ""+c1+c2;
	}
	
	
	
}
