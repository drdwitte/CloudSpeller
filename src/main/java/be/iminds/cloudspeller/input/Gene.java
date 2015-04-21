package be.iminds.cloudspeller.input;

public class Gene implements Comparable<Gene> {
	
	private String geneID;
	private String org;
	
	
	/**
	 * Constructor
	 * @param geneID
	 * @param org
	 */
	public Gene(String geneID, String org){
		this.geneID=geneID;
		this.org=org;
	}
	
	//SETTERS
	
	//GETTERS
	
	public String getID(){
		return geneID;
	}
	
	public String getOrganism(){
		return org;
	}
	
	
	
	//METHODS
	
	@Override
	public int compareTo(Gene g) {
		return this.geneID.compareTo(g.geneID);
	}
	
	public int hashCode(){
		return geneID.hashCode();
	}
	
	public boolean equals(Object obj){
		if (obj instanceof Gene) {
			Gene g = (Gene) obj;
			return geneID.equals(g.geneID);
		} else {
			return false;
		}
		
	}
	
	public String toString(){
		StringBuilder builder=new StringBuilder();
		builder.append(geneID);
		builder.append("\t");
		builder.append(org);
		return builder.toString();
	}

}