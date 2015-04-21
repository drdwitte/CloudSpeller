package be.iminds.cloudspeller.alphabets;

import java.util.HashMap;
import java.util.Map;

public class BasePairAlphabet extends Alphabet {

	public BasePairAlphabet(){
		super(generateAllChars(),generateCharMapping(),generateComplementMap());
	}
	
	private static Map<Character, Character> generateComplementMap(){
		Map<Character,Character> complementMap = new HashMap<Character, Character>();
		complementMap.put('A', 'T');
		complementMap.put('C', 'G');
		complementMap.put('G', 'C');
		complementMap.put('T', 'A');
		return complementMap;
	}

	private static Map<Character, String> generateCharMapping() {
		HashMap<Character, String> matchingChars = new HashMap<Character, String>();
		matchingChars.put('A', "A");
		matchingChars.put('C', "C");
		matchingChars.put('G', "G");
		matchingChars.put('T', "T");
		return matchingChars;
	}
	
	
	private static String generateAllChars() {
		return "ACGT";
	}

	@Override
	public CharacterIterator degenerateCharsIterator() {
		return new CharacterIterator("");
	}

	@Override
	public CharacterIterator exactCharsIterator() {
		return new CharacterIterator("ACGT");
	}

	@Override
	public int getMaxDegPerChar() {
		return 1;
	}
}
