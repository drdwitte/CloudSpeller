package be.iminds.cloudspeller.alphabets;

import java.util.Iterator;
import java.util.Map;

public abstract class Alphabet implements Iterable<Character> {
	protected final String allChars;
	protected final Map<Character, String> matchingChars;
	protected final Map<Character, Character> complementMap;

	protected Alphabet(String allChars, Map<Character, String> matchingChars
			, Map<Character, Character> complementMap) {
		this.allChars = allChars;
		this.matchingChars = matchingChars;
		this.complementMap = complementMap;
	}

	public Character generateRandomChar() {
		double random = Math.random() * allChars.length();
		int charID = (int) random;
		return allChars.charAt(charID);
	}

	@Override
	public Iterator<Character> iterator() {
		return new CharacterIterator(allChars);
	}
	
	public CharacterIterator getMatchingCharactersIterator(Character c){
		return new CharacterIterator(matchingChars.get(c));
	}
	
	public int getNumberOfMatchingCharacters(Character c){
		return matchingChars.get(c).length();
	}
	
	public boolean isDegenerate(Character c){
		return matchingChars.get(c).length()>1;
	}
	
	public CharacterIterator getAllCharsIterator(){
		return new CharacterIterator(allChars);
	}
	
	public String getAllChars() {
		return allChars;
	}

	public abstract CharacterIterator exactCharsIterator();
	public abstract CharacterIterator degenerateCharsIterator();

	public Character getComplement(char c) {
		return complementMap.get(c);
	}

	public abstract int getMaxDegPerChar();
	
}
