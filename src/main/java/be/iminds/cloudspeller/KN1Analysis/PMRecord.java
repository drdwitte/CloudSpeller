package be.iminds.cloudspeller.KN1Analysis;

import java.util.Scanner;

import be.iminds.cloudspeller.toolbox.GeneralToolbox;


public class PMRecord {

	private String motif;
	private int bls;
	private String family;
	private String gene;
	private char orientation;
	
	public PMRecord(){
		
	}
	public PMRecord(String line){
		set(line);
	}
	
	public void set(String line){	
		Scanner scanner = GeneralToolbox.generateScanner(line);
		String motifBLS	= scanner.next();
		String [] split = motifBLS.split("_");
		motif = split[0];
		bls = Integer.parseInt(split[1]);
		family		= scanner.next();
		gene		= scanner.next();
		orientation	= scanner.next().charAt(0);
				
		scanner.close();
	}
	
	
	private String concatFields(){
		return motif+bls+family+gene+orientation;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null)
			return false;
		
		if (obj instanceof PMRecord){
			PMRecord rec = (PMRecord) obj;
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
		return motif+"_"+bls+"\t"+family+"\t"+gene+"\t"+orientation;
	}
	
	/**
	 * @return the motif
	 */
	public String getMotif() {
		return motif;
	}


	/**
	 * @return the bls
	 */
	public int getBls() {
		return bls;
	}


	/**
	 * @return the family
	 */
	public String getFamily() {
		return family;
	}


	/**
	 * @return the gene
	 */
	public String getGene() {
		return gene;
	}


	/**
	 * @return the orientation
	 */
	public char getOrientation() {
		return orientation;
	}
}