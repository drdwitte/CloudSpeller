package be.iminds.cloudspeller.input;

import java.util.ArrayList;

import be.iminds.cloudspeller.alphabets.BasePairAlphabet;

public class BaseSequence extends Sequence {

	private static final BasePairAlphabet alph = new BasePairAlphabet();
	
	public BaseSequence(String s){
		super(s,alph);
	}

	public Sequence generateReverseComplement() {
		StringBuilder buffer = new StringBuilder();
		
		for (int i=0; i<super.length(); i++){
			char c = super.charAt(super.length()-1-i);
			
			Character complement = alph.getComplement(c);
			if (complement!=null){
				buffer.append(complement);
			} else {
				buffer.append(c);
			}
		}
		return new BaseSequence(buffer.toString());
	}
	
	public BaseSequence getClone(){
		return new BaseSequence(seq.toString());
	}
	
	public static ArrayList<Sequence> generateReverseComplements(ArrayList<Sequence> sequences){
		
		ArrayList<Sequence> revComps = new ArrayList<Sequence>();
		for (int i = 0; i < sequences.size(); i++) {
			Sequence currentSeq = sequences.get(i);
			BaseSequence baseSeq = (BaseSequence) currentSeq;
			revComps.add(baseSeq.generateReverseComplement());
		}
		
		return revComps;
		
	}
	
}
