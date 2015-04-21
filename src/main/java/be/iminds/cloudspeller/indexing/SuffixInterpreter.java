package be.iminds.cloudspeller.indexing;

import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.Sequence;

public class SuffixInterpreter {

	private GeneFamily gf;
	
	public SuffixInterpreter(GeneFamily gf){
		this.gf = gf;
		
	}
	
	public FullMotifMatchWithPos translateSuffixWithPos(String motif, Suffix s, int bls){
		
		FullMotifMatchWithPos match = new FullMotifMatchWithPos(translateSuffix(motif, s, bls));
				
		int numGenes = gf.getNumberOfGenes();
		Sequence seq = gf.getSequences().get(s.getSequenceID()%numGenes);
		int N = seq.length();
		int p = s.getSequencePosition();
		int k = motif.length();
		if (match.getDirection()=='+'){
			match.setMotifPosition(p-N,p-N+k-1);
		} else {
			match.setMotifPosition(-(p+1+k-1),-(p+1));
		}
		
		
		return match;
	}
	
	public static FullMotifMatchWithPos translateSuffixWithPosSingleSeq(String motif, Suffix s, String famID, String geneID, int seqLength){
		
		FullMotifMatchWithPos match = new FullMotifMatchWithPos(translateSuffixSingleSeq(motif, s, famID, geneID));

		int N = seqLength;
		int p = s.getSequencePosition();
		int k = motif.length();
		if (match.getDirection()=='+'){
			match.setMotifPosition(p-N,p-N+k-1);
		} else {
			match.setMotifPosition(-(p+1+k-1),-(p+1));
		}
	
		return match;
		
	}
	
	public static FullMotifMatch translateSuffixSingleSeq(String motif, Suffix s, String famID, String geneID){
		FullMotifMatch match = new FullMotifMatch();
		match.setGeneFamily(famID);
		match.setBls(0);
		match.setMotif(motif);
		match.setGene(geneID);
				
		char direction;
		if (s.getSequenceID() !=0){
			direction = '-';
		} else 
			direction = '+';
		match.setDirection(direction);
	
		return match;
	}
	
	
	public FullMotifMatch translateSuffix(String motif, Suffix s, int bls){
		FullMotifMatch match = new FullMotifMatch();
		
		int numGenes = gf.getNumberOfGenes();
		int seqID = s.getSequenceID()%numGenes;
				
		match.setGeneFamily(gf.getFamilyName());
		match.setBls(bls);
		
		match.setMotif(motif);
		
		String geneID = gf.getGenes().get(seqID).getID();
		match.setGene(geneID);
				
		char direction;
		if (s.getSequenceID() >= numGenes){
			direction = '-';
		} else 
			direction = '+';
		match.setDirection(direction);
		
		
		return match;
	}


}
