package be.iminds.cloudspeller.indexing;

import be.iminds.cloudspeller.input.Sequence;

import java.util.ArrayList;
import java.util.List;

public interface IndexStructureFactory {

	NodeDecorationFactory getNodeDecorationFactory();
	IndexStructure createIndexStructure(ArrayList<Sequence> seqs);
	IndexStructure createIndexStructureForSuffixes(ArrayList<Sequence> seqs, List<Suffix> suffixes);
}
