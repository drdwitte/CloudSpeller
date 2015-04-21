package be.iminds.cloudspeller.KN1Analysis;

public class RecF extends RecS {
	
	protected String geneFamily;
	protected String gene;
	
	public RecF(String geneFamily, String gene, String bindingSite, double pwmScore) {
		super(bindingSite,pwmScore);
		this.geneFamily = geneFamily;
		this.gene = gene;
	}
	
	@Override
	public String toString() {
		return geneFamily+"\t"+gene+"\t"+super.toString();
	}

	public String getGeneFamily() {
		return geneFamily;
	}

	public String getGene() {
		return gene;
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
		if (obj instanceof RecF){
			RecF record = (RecF)obj;
			return record.toString().equals(this.toString());
		}
		return false;
	}
	
	
}