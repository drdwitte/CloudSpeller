package be.iminds.cloudspeller.indexing;
import be.iminds.cloudspeller.input.Sequence;

import java.util.ArrayList;
import java.util.List;


public class GSTFactory implements IndexStructureFactory {

	private int maxDepth;
	private boolean withReverseComplements;
	private NodeDecorationFactory nodeDecoFac;
	
	public GSTFactory(int maxDepth, boolean withReverseComplements, NodeDecorationFactory fac){
		this.maxDepth=maxDepth;
		this.withReverseComplements=withReverseComplements;
		this.nodeDecoFac = fac;
		
	}
	
	@Override
	public IndexStructure createIndexStructure(ArrayList<Sequence> seqs) {
		return new GeneralizedSuffixTree(seqs, withReverseComplements,maxDepth,getNodeDecorationFactory());
	}

	@Override
	public NodeDecorationFactory getNodeDecorationFactory() {
		return nodeDecoFac;
	}

	@Override
	public IndexStructure createIndexStructureForSuffixes(
			ArrayList<Sequence> seqs, List<Suffix> suffixes) {
		return new GeneralizedSuffixTree(seqs,suffixes,maxDepth,getNodeDecorationFactory());
		
	}



}
