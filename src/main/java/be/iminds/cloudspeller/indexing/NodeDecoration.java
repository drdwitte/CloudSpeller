package be.iminds.cloudspeller.indexing;

import java.util.List;
import java.util.Set;

public interface NodeDecoration {

	void setRandomDeco();
	Set<Integer> getIDs();
	void processSuffixes(List<Suffix> suffixes);
	void joinWith(NodeDecoration nodeDecoration);
	NodeDecoration createClone();
	int toIntegerBitRepresentation();

}
