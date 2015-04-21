package be.iminds.cloudspeller.KN1Analysis;

public class RecS {

	protected String bindingSite;
	protected double pwmScore;
	
	public RecS(String bindingSite, double pwmScore) {
		this.bindingSite = bindingSite;
		this.pwmScore = pwmScore;
	}
	
	@Override
	public String toString() {
		return bindingSite+"\t"+pwmScore;
	}
	
	public String getBindingSite() {
		return bindingSite;
	}
	
	public double getPWMScore() {
		return pwmScore;
	}
	
	public void setPWMScore(double pwmScore) {
		this.pwmScore = pwmScore;
	}

	@Override
	public int hashCode() {
		return toString().hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (obj instanceof RecS){
			RecS record = (RecS)obj;
			return record.toString().equals(this.toString());
		}
		return false;
	}
	
	
}


