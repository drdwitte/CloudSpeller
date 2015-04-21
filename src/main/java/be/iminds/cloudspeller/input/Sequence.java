package be.iminds.cloudspeller.input;

import be.iminds.cloudspeller.alphabets.Alphabet;

public class Sequence {
	
	protected Alphabet sequenceAlphabet;
	protected StringBuilder seq;
	
	public Sequence(String seq, Alphabet sequenceAlphabet){
		this.seq=new StringBuilder(seq);
		this.sequenceAlphabet=sequenceAlphabet;
	}
	
	/**
	 * Copy constructor
	 * @param sequence
	 */
	public Sequence(Sequence sequence) {
		this(sequence.seq.toString(),sequence.sequenceAlphabet);
	}

	public String toString(){
		return seq.toString();
	}
	
	public boolean equals(Object obj){
		if (obj instanceof Sequence){
			Sequence s=(Sequence)obj;
			return seq.equals(s.seq);
		} else {
			return false;
		}
	}
	
	public void append(String s){
		seq.append(s);
	}
	
	public char charAt(int i){
		return seq.charAt(i);
	}
	
	public int length() {
		
		return seq.length();
	}
	
	public Alphabet getSequenceAlphabet(){
		return sequenceAlphabet;
	}

	public Sequence getClone() {
		return new Sequence(this);
	}

	public void clear() {
		seq.setLength(0);
	}
	
	public void replace(int start, int stop, String s){
		seq.replace(start,stop,s);
	}
		

}
