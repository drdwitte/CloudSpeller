package be.iminds.cloudspeller.indexing;

import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class BitSetDecoration implements NodeDecoration {

	@Override
	public String toString() {
		int [] arrayRep = new int[numberOfSequences];
		Set<Integer> sequenceSet = getIDs();
		for (Integer i : sequenceSet){
			arrayRep[i]=1;
		}
		
		StringBuilder buffer = new StringBuilder();
		for (int i=0; i<arrayRep.length; i++){
			buffer.append(arrayRep[i]);
		}
		return buffer.toString();
	}

	private int decoration;
	private int numberOfSequences;
	private int fullConservation;
	
	public BitSetDecoration(int numberOfSequences){
		this(numberOfSequences,0);
	}
	
	public BitSetDecoration(int numberOfSequences, int decoration){
		this.numberOfSequences=numberOfSequences;
		this.decoration = decoration;
		this.fullConservation = (1 << numberOfSequences) -1; //2^N-1
	}
	
	
	@Override
	public NodeDecoration createClone() {
		BitSetDecoration clone = new BitSetDecoration(this.numberOfSequences);
		clone.decoration=this.decoration;
		return clone;
	}

	@Override
	public Set<Integer> getIDs() {
		int tempValue = decoration;
		int index=0;
		Set<Integer> IDs = new HashSet<Integer>();
		while (tempValue>0){
			
			if (tempValue%2==1){
				IDs.add(index);
			}
			tempValue = tempValue >> 1;
			index++;
		}
		return IDs;
	}

	@Override
	public void joinWith(NodeDecoration nodeDecoration) {
		decoration = decoration | nodeDecoration.toIntegerBitRepresentation();
	}

	@Override
	public void processSuffixes(List<Suffix> suffixes) {
		decoration=0; //reset!
		if (suffixes.size()>numberOfSequences){ //faster version for big suffixlists
			for (Suffix s : suffixes){
				addID(s.getSequenceID());
				if (decoration == fullConservation){
					break;
				}
			}
		} else {
			for (Suffix s : suffixes){
				addID(s.getSequenceID());
				
			}
		}
	}
	
	private void addID(int id){
		int realID = id%numberOfSequences;
		decoration = decoration | (1 << realID); 
	}

	@Override
	public void setRandomDeco() {
		decoration = (int) Math.random()*(fullConservation +1);
		
	}

	@Override
	public int toIntegerBitRepresentation() {
		return decoration;
	}

}
