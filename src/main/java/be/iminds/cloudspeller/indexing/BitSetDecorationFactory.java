package be.iminds.cloudspeller.indexing;

public class BitSetDecorationFactory implements NodeDecorationFactory{

	@Override
	public NodeDecoration createNodeDecoration(int numberOfSequences) {
		return new BitSetDecoration(numberOfSequences);
	}

}