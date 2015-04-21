package be.iminds.cloudspeller.KN1Analysis;

public class ScorePair implements Comparable<ScorePair> {
	
	private String key;
	private Double value;

	public ScorePair(String key, Double value){
		this.key = key;
		this.value = value;
	}
	
	@Override
	public int compareTo(ScorePair other) {
		int comp = -this.value.compareTo(other.value);
		if (comp==0){
			return this.key.compareTo(other.key); 
		} else 
			return comp;
	}

	/**
	 * @return the key
	 */
	public String getKey() {
		return key;
	}




	/**
	 * @return the value
	 */
	public Double getValue() {
		return value;
	}




	/**
	 * @param key the key to set
	 */
	public void setKey(String key) {
		this.key = key;
	}




	/**
	 * @param value the value to set
	 */
	public void setValue(Double value) {
		this.value = value;
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
		if (!(obj instanceof ScorePair)) {
			return false;
		}
		ScorePair other = (ScorePair) obj;
		return this.toString().equals(other.toString());
	}
	
}