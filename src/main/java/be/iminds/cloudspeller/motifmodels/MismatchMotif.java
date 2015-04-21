package be.iminds.cloudspeller.motifmodels;

import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;

import org.apache.commons.lang.NotImplementedException;



public class MismatchMotif implements Motif {

	StringBuilder motif;
	int numberOfMismatches;
	
	public MismatchMotif(String m, int e) {
		this.motif=new StringBuilder(m);
		this.numberOfMismatches=e;
	}

	@Override
	public int hashCode() {
		return concatRepr().hashCode();
	}

	@Override
	public String toString() {
		return motif+"\t E="+numberOfMismatches;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof MismatchMotif){
			MismatchMotif m = (MismatchMotif) o;
			return (this.motif.equals(m.motif) && this.numberOfMismatches==m.numberOfMismatches);
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
		
	}
		
	private String concatRepr(){
		return motif.toString()+numberOfMismatches; 
	}

	

	

	@Override
	public Motif createDeepCopy() {
		return new MismatchMotif(motif.toString(),numberOfMismatches);
	}


	@Override
	public int compareTo(Motif o) {
		return this.toString().compareTo(o.toString());
	}

	@Override
	public Character charAt(int i) {
		throw new NotImplementedException();
	}

	@Override
	public MotifContent createContent() {
		throw new NotImplementedException();
	}

	@Override
	public int getGeneralizedLength() {
		throw new NotImplementedException();
	}

	@Override
	public void pop() {
		motif.deleteCharAt(motif.length()-1);		
	}

	@Override
	public Motif getComplement() {
		throw new NotImplementedException();
	}

}
