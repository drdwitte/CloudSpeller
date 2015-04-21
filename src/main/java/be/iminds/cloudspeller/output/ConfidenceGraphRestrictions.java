package be.iminds.cloudspeller.output;

public interface ConfidenceGraphRestrictions {

	public boolean checkRestrictions(ConfidenceGraph graph);

	public boolean checkRestrictions(String strValue);

	public boolean checkRestrictions(int [] F, double [] p);


}
