package be.iminds.cloudspeller.motifmodels;


import be.iminds.cloudspeller.alphabets.Alphabet;


public abstract class MotifFactory {

	public abstract Alphabet getAlphabet();
	public abstract Motif createEmptyMotif();
	
	// Motif <-> String transformation
	public abstract Motif createMotifFromString(String s);
	
	// Motif <-> BytesArray
	public abstract Motif createMotifFromBytes(int offset, byte[] total);

	public abstract byte [] createBytesRepresentation(Motif m);
	public abstract String createStringRepresentation(Motif m);
	public abstract void setMaxLength(int max);
	public abstract int getMaxLength();
	public abstract int getNumberOfBytesForMotif();
		
	public Motif createRandomMotif(int length) {
		StringBuilder randomMotif= new StringBuilder("");
		for (int i=0; i<length; i++){
			randomMotif.append(getAlphabet().generateRandomChar());
		}
		return createMotifFromString(randomMotif.toString());
	}
}
