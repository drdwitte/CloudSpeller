package be.iminds.cloudspeller.indexing;

import java.util.List;
import java.util.Set;

import org.apache.commons.lang.NotImplementedException;

public class NoDecoration implements NodeDecoration {

	@Override
	public Set<Integer> getIDs() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void processSuffixes(List<Suffix> suffixes) {
	}

	@Override
	public void setRandomDeco() {
	}

	@Override
	public void joinWith(NodeDecoration nodeDecoration) {
		throw new NotImplementedException();
	}

	@Override
	public NodeDecoration createClone() {
		throw new NotImplementedException();
	}

	@Override
	public int toIntegerBitRepresentation() {
		throw new NotImplementedException();
	}
}
