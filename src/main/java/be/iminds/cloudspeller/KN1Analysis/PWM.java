package be.iminds.cloudspeller.KN1Analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

public class PWM {

	//base order ACGT
	private double [][] matrix;
	private double [] background;
	public static double tinyProb=1e-9; //opm mismatch met tiny prob geeft ongeveer -20
	//bijdrage tot de score, gezien maximaal 
	private final int numBases=4;
	
	Map<Character,Integer> charIDMap = new HashMap<Character,Integer>();
	
	public PWM(int k){
		matrix = new double[k][];
		
		for (int i=0; i<k; i++){
			matrix[i]=new double[numBases];
		}
		
		charIDMap.put('A',0);
		charIDMap.put('C',1);
		charIDMap.put('G',2);
		charIDMap.put('T',3);
		
	}
	
	public int length(){
		return matrix.length;
	}
	
	public void setProbabilities(int index, double [] probs){
		for (int j=0; j<numBases; j++){
			if (probs[j]<tinyProb){
				matrix[index][j]=tinyProb;
			} else {
				matrix[index][j]=probs[j];
			}
		}
	}
	
	public void normalizeMatrix(){
		for (int i=0; i<matrix.length; i++){
			double sum=0.0;
			for (int j=0; j<numBases; j++){
				sum+=matrix[i][j];
			}
			for (int j=0; j<numBases; j++){
				matrix[i][j]/=sum;
			}
		}
	}
	
	public void setBackgroundModel(double [] background){
		this.background = background;
	}
	
	public double calculateMatrixMatch(String sequence){
		double score=0.0;
		for (int i=0; i<matrix.length; i++){
			
			char c = sequence.charAt(i);
			int id = charIDMap.get(c);
			
			double pi = matrix[i][id];
			double bi = background[id];
			
			//DEBUG: System.out.println(pi+"\t"+bi+"\t"+Math.log(pi/bi));
			
			score+=Math.log(pi/bi);
		}
		
		return score;
		
	}
	


	public static String createConsensus(ArrayList<String> matches, int k){
		StringBuffer consensus = new StringBuffer("");
		for (int i=0; i<k; i++){
			StringBuffer charsAtPos = new StringBuffer("");
			for (int j=0; j<matches.size(); j++){
				charsAtPos.append(matches.get(j).charAt(i));
			}
			
			char c = generateConsensusChar(charsAtPos);
			
			consensus.append(c);
		}
		
		return consensus.toString();
		
		
		
	}

	
	
	private static char generateConsensusChar(StringBuffer charsAtPos) {
		
		SortedSet<Character> chars = new TreeSet<Character>();
		for (int i=0; i<charsAtPos.length(); i++){
			chars.add(charsAtPos.charAt(i));
		}
		
		StringBuilder diffChars = new StringBuilder();
		for (Character c : chars){
			diffChars.append(c);
		}
		
		return ConsensusCharGenerator.generateConsensus(diffChars.toString());
		
		
	}
	
	

	public static void main(String [] args){
		
		String regex = "TGA..GA..GA.";
		List<StringBuilder> candSites = DegMatchesProcessor.createAllDegMatchesWithRegExp(regex,"ACGT");
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		SortedSet<MotifScore> motifScores = new TreeSet<MotifScore>();
		
		for (int i=0; i<candSites.size(); i++){
			String motif = candSites.get(i).toString();
			double score = pwm.calculateMatrixMatch(motif);
			motifScores.add(new MotifScore(motif, score));
		}
		
		for (MotifScore m : motifScores){
			System.out.println(m);
		}
		
		
		
		
		
		String worst = "AACAAACACACA";
		System.out.println("Worst motif: "+worst + "\t" + pwm.calculateMatrixMatch(worst));
		
		
		
		System.out.println("Now lets investigate consensus motifs");
		
		ArrayList<String> matches = new ArrayList<String>();
		for (MotifScore m : motifScores){
			matches.add(m.getMotif());
			String consensus = createConsensus(matches,12);
			
			System.out.println(new MotifScore(consensus,m.getScore()));
			
		}
		
		
		
		
		
		
		System.out.println("done");
	}

		
	
}

class MotifScore implements Comparable<MotifScore> {
	/**
	 * @return the motif
	 */
	public String getMotif() {
		return motif;
	}

	/**
	 * @param motif the motif to set
	 */
	public void setMotif(String motif) {
		this.motif = motif;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	private String motif;
	private double score;
	
	public MotifScore(String motif, double score){
		this.motif=motif;
		this.score=score;
	}

	@Override
	public int compareTo(MotifScore o) {
		
		if (this.score<o.score){
			return +1;
		} else if (this.score>o.score){
			return -1;
		} else {
			return this.motif.compareTo(o.motif); 
		}
		
	}
	
	public String toString(){
		return motif + "\t" + score;
	}
}

class ConsensusCharGenerator {

	/*
	matchingChars.put('M', "AC");
	matchingChars.put('R', "AG");
	matchingChars.put('W', "AT");
	matchingChars.put('S', "CG");
	matchingChars.put('Y', "CT");
	matchingChars.put('K', "GT");
	*/
	
	public static char generateConsensus(String chars) {
		
		if (chars.length()==1){
			return chars.charAt(0);
		} else if (chars.length()==2){
			if (chars.charAt(0)=='A'){  
				if (chars.charAt(1) == 'C'){ //AC
					return 'M';
				} else if (chars.charAt(1) == 'G'){ //AG 
					return 'R';
				} else { //AT
					return 'W';
				}
			} else if (chars.charAt(0)=='C'){ 
				if (chars.charAt(1) == 'G'){ //CG 
					return 'S';
				} else { //CT
					return 'Y';
				}
			} else { //GT
				return 'K';
			
			} 
		}
		return 'N';
		
	}
	
}

