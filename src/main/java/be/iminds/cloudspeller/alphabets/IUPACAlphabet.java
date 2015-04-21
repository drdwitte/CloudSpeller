package be.iminds.cloudspeller.alphabets;

import java.util.HashMap;
import java.util.Map;

public class IUPACAlphabet extends Alphabet {
	
	public String degChars;
	private int maxDegPerChar=4;

	public IUPACAlphabet(IUPACType type) {
		super(generateAllChars(type), generateCharMapping(type),generateComplementMap(type));
		
		if (type.equals(IUPACType.BASEPAIRS)){
			maxDegPerChar=1;
		}
		
		degChars = allChars.substring(4);

	}
	
	private static Map<Character, Character> generateComplementMap(IUPACType type){
	
		Map<Character,Character> complementMap = new HashMap<Character, Character>();

		switch (type) { // NOTE switch without breaks!
		case FULL:
			complementMap.put('B','V');
			complementMap.put('D','H');
			complementMap.put('H','D');
			complementMap.put('V','B');
		
		case TWOFOLDSANDN:
			complementMap.put('M','K');
			complementMap.put('R','Y');
			complementMap.put('W','W'); 
			complementMap.put('S','S');
			complementMap.put('Y','R');
			complementMap.put('K','M');
			
		case DONTCARES:
			complementMap.put('N','N');
			
		case BASEPAIRS:
			complementMap.put('A','T');
			complementMap.put('C','G');
			complementMap.put('G','C');
			complementMap.put('T','A');
		}
		return complementMap;
	}


	private static Map<Character, String> generateCharMapping(IUPACType type) {
		HashMap<Character, String> matchingChars = new HashMap<Character, String>();
		
		switch (type) { // NOTE switch without breaks!
		case FULL:
			matchingChars.put('B', "CGT");
			matchingChars.put('D', "AGT");
			matchingChars.put('H', "ACT");
			matchingChars.put('V', "ACG");
			
		case TWOFOLDSANDN:
			matchingChars.put('M', "AC");
			matchingChars.put('R', "AG");
			matchingChars.put('W', "AT");
			matchingChars.put('S', "CG");
			matchingChars.put('Y', "CT");
			matchingChars.put('K', "GT");
			
		case DONTCARES:
			matchingChars.put('N', "ACGT");
			
		case BASEPAIRS:
			matchingChars.put('A', "A");
			matchingChars.put('C', "C");
			matchingChars.put('G', "G");
			matchingChars.put('T', "T");
		}
		return matchingChars;
	}

	private static String generateAllChars(IUPACType type) {
		StringBuilder alphabetCharacters = new StringBuilder();
		
		switch (type) { // NOTE switch without breaks!
		case FULL:
			alphabetCharacters.insert(0, "BDHV"); // !A!C!G!T
		case TWOFOLDSANDN:
			alphabetCharacters.insert(0, "MRWSYK"); // AC AG AT CG CT GT
		case DONTCARES:
			alphabetCharacters.insert(0, "N");
		case BASEPAIRS:
			alphabetCharacters.insert(0, "ACGT");
		}
		
		return alphabetCharacters.toString();
	}

	public static enum IUPACType {
		BASEPAIRS, DONTCARES, TWOFOLDSANDN, FULL
	}



	@Override
	public CharacterIterator degenerateCharsIterator() {
		return new CharacterIterator(degChars);
	}
	
	

	@Override
	public CharacterIterator exactCharsIterator() {
		return new CharacterIterator("ACGT");
	}
	
	public String getDegChars(){
		return degChars;
	}

	@Override
	public int getMaxDegPerChar() {
		return maxDegPerChar;
	}
}
