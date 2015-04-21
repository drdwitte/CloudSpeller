package be.iminds.cloudspeller.motifmodels;

import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;

import org.apache.commons.lang.NotImplementedException;




public class BasePairMotif implements Motif {

	private StringBuilder motif;

	/**
	 * Constructor
	 * @param string
	 */
	public BasePairMotif(String string) {
		this.motif=new StringBuilder(string);
	}
	
	private BasePairMotif(){}
	
	//GETTERS
	@Override
	public int length() {
		return motif.length();
	}
	
	//SETTERS
	
	//METHODS

	@Override
	public boolean equals(Object o) {
		if (o instanceof BasePairMotif){
			BasePairMotif bp= (BasePairMotif) o;
			return motif.toString().equals(bp.motif.toString());
		} else 
			return false;
	}
	
	@Override
	public int hashCode() {
		return motif.toString().hashCode();
	}
		
	@Override
	public String toString(){
		return motif.toString();
	}

	@Override
	public void append(Character c) {
		motif.append(c);
	}

	@Override
	public Motif createDeepCopy() {
		BasePairMotif bp = new BasePairMotif();
		bp.motif=new StringBuilder(this.motif);
		return bp;
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
