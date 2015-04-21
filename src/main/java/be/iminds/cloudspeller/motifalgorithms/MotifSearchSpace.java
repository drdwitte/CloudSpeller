package be.iminds.cloudspeller.motifalgorithms;

import be.iminds.cloudspeller.alphabets.Alphabet;

public class MotifSearchSpace {
	
	private int minLength;
	private int maxLength;
	private int maxDegeneratePositions;
	private Alphabet alphabet;

	public MotifSearchSpace(int kmin, int kmax, int maxDegPos, Alphabet alph) {
		this.minLength=kmin;
		this.maxLength=kmax;
		this.maxDegeneratePositions=maxDegPos;
		this.alphabet=alph;
		
	}

	public int getMinLength() {
		return minLength;
	}

	public int getMaxLength() {
		return maxLength;
	}

	public int getMaxNumberOfDegeneratePositions() {
		return maxDegeneratePositions;
	}

	public Alphabet getAlphabet() {
		return alphabet;
	}
	
	public String toString(){
		StringBuilder b = new StringBuilder();
		b.append("(");
		b.append(minLength);
		b.append(",");
		b.append(maxLength);
		b.append(",");
		b.append(maxDegeneratePositions);
		b.append(")");
		
		return b.toString();
	}

	public int getMaximumDegeneracy() {
		int maxDeg = 1;
		for (int i=0; i<maxDegeneratePositions; i++){
			maxDeg*=alphabet.getMaxDegPerChar();
		}
		return maxDeg;
		
	}

}
