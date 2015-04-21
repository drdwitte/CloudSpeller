package be.iminds.cloudspeller.indexing;

import java.util.ArrayList;
import java.util.List;

import be.iminds.cloudspeller.motifmodels.Motif;

public interface ISMonkey {

	public Motif getMotifTrail();
	public boolean hasMatches();
	public List<Suffix> grabSuffixes();
	public ArrayList<NodeDecoration> grabInternalNodeInfo();
	public ISMonkey createClone();
	public void jumpTo(Character c);
	public void backtrack();
	//public ISMonkey createDegenerateMonkey(Character degenerateExtension,
	//		List<ISMonkey> exactMonkeys);
	public NodeDecoration getJointNodeInfo();
}
