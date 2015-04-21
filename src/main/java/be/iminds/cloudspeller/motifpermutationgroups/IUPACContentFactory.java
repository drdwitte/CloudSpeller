package be.iminds.cloudspeller.motifpermutationgroups;

import org.apache.hadoop.io.BytesWritable;

import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;


public class IUPACContentFactory implements MotifContentFactory {
	
	private MotifFactory motifFactory; 

	public IUPACContentFactory(MotifFactory motifFactory){
		this.motifFactory=motifFactory;
	}
	
	@Override
	public MotifContent createMotifContentFromBytes(BytesWritable b) {
		Motif m = motifFactory.createMotifFromBytes(0,b.getBytes());
		return new IUPACContent(m.toString()); 
								//use string since already sorted!
	}

	@Override
	public MotifContent createMotifContentFromString(String s) {
		return new IUPACContent(s);
	}

	@Override
	public MotifContent createContent(Motif m) {
		return new IUPACContent(m);
	}

	@Override
	public BytesWritable createBytesRepresentation(MotifContent groupID) {
		Motif m = motifFactory.createMotifFromString(groupID.toString());
		return new BytesWritable(motifFactory.createBytesRepresentation(m));
	}

	@Override
	public String createStringRepresentation(MotifContent groupID) {
		return groupID.toString();
	}
	
}