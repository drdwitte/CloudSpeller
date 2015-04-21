package be.iminds.cloudspeller.alphabets;

import java.util.Iterator;

public class CharacterIterator implements Iterator<Character> {

	private int pos = 0;
	private String text;

	public CharacterIterator(String text) {
		this.text = text;
	}

	@Override
	public boolean hasNext() {
		return pos < text.length();
	}

	@Override
	public Character next() {
		return text.charAt(pos++);
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException(
				"Removing alphabet chars not allowed!");

	}
	
	public int getAlphabetSize(){
		return text.length();
	}
}
