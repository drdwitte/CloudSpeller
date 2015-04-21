package be.iminds.cloudspeller.indexing;

public class NoDecorationFactory implements NodeDecorationFactory {

	@Override
	public NodeDecoration createNodeDecoration(int numberOfSequences) {
		return new NoDecoration();
	}

}
