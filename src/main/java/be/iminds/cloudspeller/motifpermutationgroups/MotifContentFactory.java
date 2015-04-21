package be.iminds.cloudspeller.motifpermutationgroups;

import org.apache.hadoop.io.BytesWritable;

import be.iminds.cloudspeller.motifmodels.Motif;



public interface MotifContentFactory {

	
	
	//streaming functionalities
	public MotifContent createMotifContentFromString(String s);
	public MotifContent createMotifContentFromBytes(BytesWritable b);
	
	public MotifContent createContent(Motif m);
	public String createStringRepresentation(MotifContent groupID);
	public BytesWritable createBytesRepresentation(MotifContent groupID);
	
}
