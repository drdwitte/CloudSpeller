package be.iminds.cloudspeller.KN1Analysis;

import java.util.Locale;
import java.util.Scanner;

public class DeNovoRecord {
	

	private String motif;
	private int blsMin;
	private int numFamilies;
	private double confidence;
	private double pwmScore;
	
	
	public DeNovoRecord(String line){
		Scanner scanner = new Scanner(line);
		scanner.useLocale(Locale.US);
		motif		= scanner.next();
		blsMin		= scanner.nextInt();
		numFamilies	= scanner.nextInt();
		confidence	= scanner.nextDouble();
		pwmScore	= scanner.nextDouble();
		scanner.close();
	}
	
	
	private String concatFields(){
		return motif+blsMin+numFamilies+confidence+pwmScore;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null)
			return false;
		
		if (obj instanceof DeNovoRecord){
			DeNovoRecord rec = (DeNovoRecord) obj;
			return this.concatFields().equals(rec.concatFields());
		}
		return false;
	}
	
	@Override
	public int hashCode() {
		return this.concatFields().hashCode();
	}
	
	@Override 
	public String toString() {
		return motif+"\t"+blsMin+"\t"+numFamilies+"\t"+confidence+"\t"+pwmScore;
	}
	
	/**
	 * @return the motif
	 */
	public String getMotif() {
		return motif;
	}


	/**
	 * @return the blsMin
	 */
	public int getBlsMin() {
		return blsMin;
	}


	/**
	 * @return the numFamilies
	 */
	public int getNumFamilies() {
		return numFamilies;
	}


	/**
	 * @return the confidence
	 */
	public double getConfidence() {
		return confidence;
	}


	/**
	 * @return the pwmScore
	 */
	public double getPwmScore() {
		return pwmScore;
	}
}