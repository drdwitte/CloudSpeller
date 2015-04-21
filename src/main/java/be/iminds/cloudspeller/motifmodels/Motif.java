package be.iminds.cloudspeller.motifmodels;

import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;

public interface Motif extends Comparable<Motif> {

	public int hashCode();
	public boolean equals(Object o);
	public String toString();
	public int length();
	public int getGeneralizedLength();
	public void append(Character c);
	public Motif createDeepCopy();
	public Character charAt(int i);
	public MotifContent createContent();
	public void pop();
	public Motif getComplement();
}
