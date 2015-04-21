package be.iminds.cloudspeller.postprocessing;

import java.util.HashSet;
import java.util.Set;


import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.CharacterIterator;

/**
 * character distance calculator based on jacard distance, note that d(char,indel)=BIG_DISTANCE
 * @author ddewitte
 *
 */

public class JacardDistanceCalculator implements CharacterDistanceCalculator {

	public static final char indelCharBulk = '-';
	public static final char indelCharBoundary = '.';
	protected static final double BIG_DISTANCE = 1000;
	protected Alphabet alphabet;
			
	/*
	 * Constructor
	 */
	public JacardDistanceCalculator(Alphabet alphabet){
		setAlphabet(alphabet);
	}
	
	//GETTERS 
	@Override
	public char getBoundaryIndelCharacter() {
		return indelCharBoundary;
	}

	@Override
	public char getBulkIndelCharacter() {
		return indelCharBulk;
	}
	
	//SETTERS
	public void setAlphabet(Alphabet alphabet){
		this.alphabet = alphabet;
	}
	
	
	
	//METHODS
	
	@Override
	public double calculateDistance(char c1, char c2){
		
		if (c1 == indelCharBulk || c2 == indelCharBulk){
			return BIG_DISTANCE;
		}
		
		int intersection;
		int union;
		
		if (c1 == indelCharBoundary || c2 == indelCharBoundary) {
			
			if (c1 == indelCharBoundary){
				
				if (c2 == indelCharBoundary){
					return 0.0;
				} else {
					intersection = alphabet.getNumberOfMatchingCharacters(c2);
				}
			} else {
				intersection = alphabet.getNumberOfMatchingCharacters(c1);
			}
			union = alphabet.exactCharsIterator().getAlphabetSize();
					
			
		} else {
			intersection = calcIntersection(c1,c2);
			if (intersection==0){
				return 1.0;
			}
			union = calcUnion(c1, c2);
		}
		
		return 1.0 - 1.0*intersection/union;
	}
	
	protected int calcUnion(char c1, char c2) {
		
		CharacterIterator it1 = alphabet.getMatchingCharactersIterator(c1);
		CharacterIterator it2 = alphabet.getMatchingCharactersIterator(c2);
		
		Set<Character> uniqueChars = new HashSet<Character>();
		
		while (it1.hasNext()){
			uniqueChars.add(it1.next());
		}
		
		while (it2.hasNext()){
			uniqueChars.add(it2.next());
		}
		
		return uniqueChars.size();
		
	}
	
	protected int calcIntersection(char c1, char c2) {
		CharacterIterator it1 = alphabet.getMatchingCharactersIterator(c1);
		int numberOfMatches=0;
		while (it1.hasNext()){
			char alphChar1=it1.next();
			CharacterIterator it2 = alphabet.getMatchingCharactersIterator(c2);
				while (it2.hasNext()){
				char alphChar2=it2.next();
				
				if (alphChar1==alphChar2){
					numberOfMatches++;
				}
			}
		}
		
		return numberOfMatches;
	}

	
}
