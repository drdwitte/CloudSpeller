package be.iminds.cloudspeller.motifmodels;



import be.iminds.cloudspeller.motifpermutationgroups.IUPACContent;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;
import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.CharacterIterator;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

public class IUPACMotif implements Motif {

	private static final Alphabet alphabet = new IUPACAlphabet(IUPACType.FULL);
	private StringBuilder motif;
	private int numberOfDegeneratePositions;

	private IUPACMotif(){
		
	}
	
	public IUPACMotif(String s){
		this.motif=new StringBuilder(s);
		calculateNumberOfDegeneratePositions(s);
	}
	
	public IUPACMotif(String s, int numberOfDegPos){
		this.motif = new StringBuilder(s);
		this.numberOfDegeneratePositions=numberOfDegPos;
	}
	
	public int calculateDegeneracy(String s){
		CharacterIterator iterator = new CharacterIterator(s);
		int degeneracy = 1;
		while (iterator.hasNext()){
			degeneracy*=alphabet.getNumberOfMatchingCharacters(iterator.next());
		}
		return degeneracy;
	}
	
	public void calculateNumberOfDegeneratePositions(String s){
		CharacterIterator iterator = new CharacterIterator(s);
		numberOfDegeneratePositions = 0;
		while (iterator.hasNext()){
			numberOfDegeneratePositions+=addDegPositionContribution(iterator.next());
		}
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof IUPACMotif){
			IUPACMotif m =(IUPACMotif) o;
			return m.motif.toString().equals(motif.toString());
		}
		return false;
	}

	@Override
	public int length() {
		return motif.length();
	}

	@Override
	public void append(Character c) {
		motif.append(c);
		numberOfDegeneratePositions+=addDegPositionContribution(c);
		
	}
	
	
	private int addDegPositionContribution(Character c){
		return (alphabet.getNumberOfMatchingCharacters(c)>1)?1:0;
	}

	@Override
	public int hashCode() {
		return motif.toString().hashCode();
	}

	@Override
	public String toString() {
		return motif.toString();
	}

	@Override
	public Character charAt(int i) {
		return motif.charAt(i);
	}

	@Override
	public Motif createDeepCopy() {
		IUPACMotif m = new IUPACMotif();
		m.motif=new StringBuilder(this.motif);
		m.numberOfDegeneratePositions=this.numberOfDegeneratePositions;
		return m;
	}

	@Override
	public int compareTo(Motif o) {
		return this.toString().compareTo(o.toString());
	}

	public int numberOfDegPositions() {
		return numberOfDegeneratePositions;
	}

	@Override
	public MotifContent createContent() {
		return new IUPACContent(this);
	}

	@Override
	public int getGeneralizedLength() {
		return this.length()-this.numberOfDegeneratePositions/2;
	}

	@Override
	public void pop() {
		numberOfDegeneratePositions -= addDegPositionContribution(motif.charAt(motif.length()-1));
		motif.deleteCharAt(motif.length()-1);
		
	}

	@Override
	public Motif getComplement() {
		StringBuilder s = new StringBuilder();
		for (int i=motif.length()-1; i>=0; i--){
			s.append(alphabet.getComplement(motif.charAt(i)));
		}
		return new IUPACMotif(s.toString(),numberOfDegeneratePositions);
	}

}
