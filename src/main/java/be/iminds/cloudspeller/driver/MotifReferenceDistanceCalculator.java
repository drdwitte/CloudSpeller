package be.iminds.cloudspeller.driver;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.postprocessing.JacardDistanceCalculator;

public class MotifReferenceDistanceCalculator extends JacardDistanceCalculator {

	/*
	 * Constructor
	 */
	public MotifReferenceDistanceCalculator(Alphabet alphabet) {
		super(alphabet);
	}
	
	
	/**
	 * d(A|C,N)=1, d(A,A|C)=1, d(A,N)=2, d(A|C,A|T)=1.5
	 */
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
					System.out.println(0.0);
					return 0.0;
				} else {
					intersection = alphabet.getNumberOfMatchingCharacters(c2);
				}
			} else {
				System.out.println(alphabet.getAllChars());
				intersection = alphabet.getNumberOfMatchingCharacters(c1);
			}
			union = alphabet.exactCharsIterator().getAlphabetSize();
					
			
		} else {
			intersection = calcIntersection(c1,c2);
			if (intersection==0){
				return BIG_DISTANCE;
			}
			union = calcUnion(c1, c2);
		}
		
		
		return -0.5+0.5*union/intersection;
	}




}
