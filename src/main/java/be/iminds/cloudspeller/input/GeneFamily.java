package be.iminds.cloudspeller.input;

import java.util.*;
import java.io.*;

public class GeneFamily {

	//public static final char GFSTART = '>';
	//public static boolean addReverseComplements=false;
	private static String generalNewickString;
	private String famID;
	private StringBuilder newick;
	private SortedMap<Gene,BaseSequence> geneSeq=new TreeMap<Gene,BaseSequence>();
		
	
	/**
	 * Gene Family constructor (only assigns name)
	 * @param famID Gene family identifier
	 */
	public GeneFamily(String famID){
		if (generalNewickString!=null){
			
		}
		this.famID=famID;
	}
	
	public boolean isIntialized(){
		return geneSeq.size()>0 && famID!=null && newick!=null;
	}
	
	/**
	 * Gene Family constructor from buffered reader
	 * @param in Buffered reader to extract gene family information from
	 * @throws IOException
	 */
	public GeneFamily(BufferedReader in) throws IOException {
		String line;
		while ((line=in.readLine())!=null){
			this.famID=line; 
			this.newick=new StringBuilder(in.readLine());
			int numGenes=Integer.parseInt(in.readLine());
			for (int i=0; i<numGenes; i++){
				line=in.readLine();
				Scanner scan=new Scanner(line);
				String geneID=scan.next();
				String orgName=scan.next();
				Gene g=new Gene(geneID,orgName);
				String seqStr=in.readLine();
				BaseSequence seq=new BaseSequence(seqStr);
				geneSeq.put(g, seq);
				
				scan.close();
			}
			break;
		}
	}
	

	//SETTERS
	/**
	 * Set newick string externally
	 * @param newick newick string to set
	 */
	public void setNewick(String newick){
		this.newick=new StringBuilder(newick);
	}
	
	/**
	 * General newick string containing the species names, which can be overridden with gene names
	 * @param genNewick
	 */
	public static void setGeneralNewick(String genNewick){
		generalNewickString=genNewick;
	}
	
	
	//GETTERS
	public String getFamilyName(){
		return famID;
	}
	
	public int getNumberOfGenes(){
		return geneSeq.size();
	}
	
	public ArrayList<Sequence> getSequences() {
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();
		for (Map.Entry<Gene,BaseSequence> entry : geneSeq.entrySet()){
			sequences.add(entry.getValue());
		}
		return sequences;
	}
	
	public ArrayList<Gene> getGenes() {
		ArrayList<Gene> sequences = new ArrayList<Gene>();
		for (Map.Entry<Gene,BaseSequence> entry : geneSeq.entrySet()){
			sequences.add(entry.getKey());
		}
		return sequences;
	}
	
	public StringBuilder getNewick(){
		return newick;
	}
	
	
	
	
	//METHODS
	
	/**
	 * @return true if both the sequences and the newickstring have been initialized
	 */
	public boolean isInitialized(){
		return (geneSeq.size()!=0 && newick.length()!=0);
	}
	
	public int hashCode(){
		return famID.hashCode();
	}
	

	/**
	 * Generates unique newick string which is not necessarily normalized
	 * @param paralogBranchLength
	 */
	public void generateNewick(double paralogBranchLength){
		if (paralogBranchLength<=0.0){
			throw new UnsupportedOperationException( "Paralog branch lenght must be strictly > 0!");
		}
		this.newick=new StringBuilder(generalNewickString);
		
		//how many genes per org?
		Map<String,HashSet<String> > orgMap= new HashMap<String,HashSet<String>>();
		
		for(Map.Entry<Gene,BaseSequence> entry : geneSeq.entrySet()){
			Gene g=entry.getKey();
			HashSet<String> gs=orgMap.get(g.getOrganism());
			if (gs!=null){
				gs.add(g.getID());
			} else {
				HashSet<String> newGs=new HashSet<String>();
				newGs.add(g.getID());
				orgMap.put(g.getOrganism(),newGs);
			}
		}
		
		for (Map.Entry<String,HashSet<String>> entry : orgMap.entrySet()){
			HashSet<String> genesFromOneOrg=entry.getValue();
			
			if (genesFromOneOrg.size()<=1){
				if (genesFromOneOrg.size()==0){
					System.out.println("Not all genes present in: "+famID);
					newick.setLength(0);
				}
				else {
					String org=entry.getKey();
					int start=newick.indexOf(org);
					if (start<0){
						System.out.println("Newick error: org not found!");
					}
					int stop=start+org.length();
					newick.replace(start,stop, genesFromOneOrg.iterator().next());
				}
			} else {
				String org=entry.getKey();
				int start=newick.indexOf(org);
				if (start<0){
					System.out.println("Newick error: org not found!");
				}
				int stop=start+org.length();
				
				StringBuilder temp= new StringBuilder();
				temp.append("(");
				
				String paralogBranch=":"+paralogBranchLength;
				
				Iterator<String> itG=genesFromOneOrg.iterator();
				while (itG.hasNext()){
					temp.append(itG.next());
					temp.append(paralogBranch);
					temp.append(",");
				}
				
				temp.replace(temp.length()-1,temp.length(), ")");
				newick.replace(start,stop, temp.toString());
				
			}
			
		}
		
	}
	
	
	public void addGeneSeq(Gene g, BaseSequence seq){
		geneSeq.put(g,seq);
	}
	
	public String toString(){
		StringBuilder famStringB=new StringBuilder();
		
		//famStringB.append(""+GFSTART);
		famStringB.append(famID);
		famStringB.append("\n");
		famStringB.append(newick);
		famStringB.append("\n");
		famStringB.append(geneSeq.size());
		famStringB.append("\n");
		
		for (Map.Entry<Gene,BaseSequence> entry : geneSeq.entrySet()){
	
			Gene g=entry.getKey();
			famStringB.append(g);
			famStringB.append("\n");
			famStringB.append(entry.getValue().toString());
			famStringB.append("\n");
		}

		return famStringB.toString();
	}

	public Sequence getSequence(Gene g) {

		return geneSeq.get(g);
	}

	
		
	




	
	
	
}