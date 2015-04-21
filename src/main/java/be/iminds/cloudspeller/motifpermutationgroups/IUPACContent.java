package be.iminds.cloudspeller.motifpermutationgroups;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;



import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;


public class IUPACContent implements MotifContent {

	private String contentString;
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof IUPACContent){
			return ((IUPACContent) obj).contentString.equals(this.contentString);
		}
		return false;
	}
	
	@Override
	public int hashCode() {
		return contentString.hashCode();
	}
	
	@Override
	public String toString() {
		return contentString;
	}
	
	public IUPACContent(Motif m){
	
		setContentFromStr(m.toString());
	}
	
	public void setContentFromStr(String motifStr){
		char[] chars = motifStr.toCharArray();
	    Arrays.sort(chars);
	    this.contentString = new String(chars);
	}
	
	public IUPACContent(String sorted){
		this.contentString = sorted;
	}

	@Override
	/**
	 * Convert contentstring into a List<Character> and use collections.shuffle to generate
	 * random shuffles = permutation group
	 * @param numberOfTries Number of shuffles (not that Set size is not determined beforehand) 
	 */
	public Set<Motif> createPermutationGroup(int numberOfTries) {
		MotifFactory factory= new IUPACFactory(IUPACType.FULL);
		
		Set<Motif> backgroundMotifs = new HashSet<Motif>();
		List<Character> chars=charListFromString(contentString);

		for (int i=0; i<numberOfTries; i++){
			Collections.shuffle(chars);
			Motif motif=factory.createEmptyMotif();
			Iterator<Character> it=chars.iterator();
			while (it.hasNext()){
				motif.append(it.next());
			}
			backgroundMotifs.add(motif);
		}
		return backgroundMotifs;
	}
	
	public static List<Character> charListFromString(String s){
		List<Character> list = new ArrayList<Character>();
		for (int i=0; i<s.length(); i++){
			list.add(s.charAt(i));
		}
		return list;
	}

}




