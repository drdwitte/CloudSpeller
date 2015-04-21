package be.iminds.cloudspeller.indexing;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SequenceIDSet implements NodeDecoration {

	private int numberOfSequences;
	private Set<Integer> sequenceSet=new HashSet<Integer>();
	
	/**
	 * Constructor
	 * @param numberOfSequences
	 */
	public SequenceIDSet(int numberOfSequences) {
		this.numberOfSequences=numberOfSequences;
	}

	//SETTERS
	
	@Override
	public void setRandomDeco() {
		sequenceSet.clear();
		int numSeq=1+(int)(Math.random()*(numberOfSequences));
		
		for (int i=0; i<numSeq; i++){
			int randSeq=(int)(Math.random()*numberOfSequences);
			sequenceSet.add(randSeq);
		}
	}

	@Override
	public Set<Integer> getIDs() {
		return sequenceSet;
	}

	@Override
	public void processSuffixes(List<Suffix> suffixes) {
		sequenceSet.clear();
		if (suffixes.size()>numberOfSequences){ //faster version for big suffixlists
			for (Suffix s : suffixes){
				addID(s.getSequenceID());
				if (sequenceSet.size()==numberOfSequences){
					break;
				}
			}
		} else {
			for (Suffix s : suffixes){
				addID(s.getSequenceID());
				
			}
		}
	
	}



	private void addID(int sequenceID) {
		sequenceSet.add(sequenceID%numberOfSequences);
	}

	@Override
	public void joinWith(NodeDecoration nodeDecoration) {
		
		for (Integer i : nodeDecoration.getIDs()){
			sequenceSet.add(i);
		}
	}

	@Override
	public String toString() {
		int [] arrayRep = new int[numberOfSequences];
		for (Integer i : sequenceSet){
			arrayRep[i]=1;
		}
		
		StringBuilder buffer = new StringBuilder();
		for (int i=0; i<arrayRep.length; i++){
			buffer.append(arrayRep[i]);
		}
		return buffer.toString();
		
	}

	@Override
	public NodeDecoration createClone() {
		NodeDecoration clone = new SequenceIDSet(this.numberOfSequences);
		clone.joinWith(this);
		return clone;
	}

	@Override
	public int toIntegerBitRepresentation() {
		int repr=0;
		
		for (Integer i : sequenceSet){
			repr+= 1 << i;
		}
		return repr;
	}
}


