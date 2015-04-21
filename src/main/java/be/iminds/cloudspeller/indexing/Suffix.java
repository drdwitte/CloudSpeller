package be.iminds.cloudspeller.indexing;

public class Suffix {
	
	private int sequenceID;
	private int seqPos;
	
	public Suffix(int seqID, int seqP){
		setSequenceID(seqID);
		setSequencePosition(seqP);
	}

	public int getSequenceID() {
		return sequenceID;
	}

	public void setSequenceID(int sequenceID) {
		this.sequenceID = sequenceID;
	}

	public int getSequencePosition() {
		return seqPos;
	}

	public void setSequencePosition(int seqPos) {
		this.seqPos = seqPos;
	}

	@Override
	public String toString() {
		return "("+sequenceID+","+seqPos+")";
	}

}
