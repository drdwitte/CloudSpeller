package be.iminds.cloudspeller.motifmodels;

import org.apache.commons.lang.NotImplementedException;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.BasePairAlphabet;

public class MMFactory extends MotifFactory {

	private Alphabet alphabet = new BasePairAlphabet();
	
	@Override
	public Motif createEmptyMotif() {
		throw new NotImplementedException();
	}

	@Override
	public Motif createMotifFromBytes(int offset, byte [] b) {
		throw new NotImplementedException();
	}

	@Override
	public Motif createMotifFromString(String s) {
		throw new NotImplementedException();
	}

	@Override
	public Alphabet getAlphabet() {
		return alphabet;
	}

	@Override
	public int getNumberOfBytesForMotif() {
		throw new NotImplementedException();
	}

	@Override
	public byte[] createBytesRepresentation(Motif m) {
		throw new NotImplementedException();
	}

	@Override
	public String createStringRepresentation(Motif m) {
		throw new NotImplementedException();
	}

	@Override
	public int getMaxLength() {
		throw new NotImplementedException();
	}

	@Override
	public void setMaxLength(int max) {
		throw new NotImplementedException();
	}
}
