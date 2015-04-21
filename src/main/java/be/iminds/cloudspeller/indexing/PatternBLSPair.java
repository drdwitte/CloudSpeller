package be.iminds.cloudspeller.indexing;

import java.util.Scanner;

public class PatternBLSPair {

	private String pattern;
	private int bls;
	
	public PatternBLSPair(String line){
		Scanner scan = new Scanner(line);
		setPattern(scan.next());
		setBls(scan.nextInt());
		scan.close();
	}
	
	public PatternBLSPair(String pattern, int bls){
		setPattern(pattern);
		setBls(bls);
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PatternBLSPair){
			return ((PatternBLSPair)obj).pattern.equals(this.pattern)
					&& ((PatternBLSPair)obj).bls==this.bls;
		} 
		return false;
	}

	@Override
	public int hashCode() {
		return toString().hashCode();
	}

	@Override
	public String toString() {
		return pattern+"_"+bls;
	}

	public String getPattern() {
		return pattern;
	}

	public void setPattern(String pattern) {
		this.pattern = pattern;
	}

	public int getBls() {
		return bls;
	}

	public void setBls(int bls) {
		this.bls = bls;
	}
	
	
	
}
