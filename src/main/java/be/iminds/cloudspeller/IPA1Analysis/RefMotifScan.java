package be.iminds.cloudspeller.IPA1Analysis;

import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.Sequence;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.KN1Analysis.RecS;



public class RefMotifScan {
	
	public static Map<Gene, Set<RecS>> getAllGenesWithRefMotif(
			Map<Gene, Sequence> genes, String referenceMotif) {
		int c=0;
		Map<Gene, Set<RecS>> recordMap = new HashMap<Gene, Set<RecS>>();
		
		for (Map.Entry<Gene,Sequence> e : genes.entrySet()){
			c++;
			if (c%100==0){
				System.out.println("promoters scanned: "+c);
			}
			Sequence seq = e.getValue();
			List<Suffix> matches = IPA1Toolbox.giveMatchesOfPatternInSequence(referenceMotif, seq);
						
			if (matches!=null){
				Set<RecS> records = new HashSet<RecS>();
				for (Suffix m : matches){
					String bs = IPA1Toolbox.getActualBS(m, seq);
					records.add(new RecS(bs,-1));
				}
				
				recordMap.put(e.getKey(),records);
			}
		
		}
		return recordMap;
		
	}

	public static Map<Gene, Sequence> getAllGeneDataFromMonocots(String species) {
		
		Map<Gene, Sequence> geneSeqs = new HashMap<Gene, Sequence>();
		
		for (Map.Entry<String,GeneFamily> e : IPA1Toolbox.getFamilies().entrySet()){
			
			ArrayList<Gene> genes = e.getValue().getGenes();
			ArrayList<Sequence> seqs = e.getValue().getSequences();
			
			for (int i=0; i<genes.size(); i++){
				Gene g = genes.get(i);
				if (g.getOrganism().equals(species)){ //is species gene?
					
					geneSeqs.put(genes.get(i),seqs.get(i));
				}
			}
		} 
		return geneSeqs;
	}
	
	public static void main (String [] args) throws IOException{
		
		IPA1Toolbox.setDatasetFromFile();
			
		
		Map<Gene,Sequence> genes = getAllGeneDataFromMonocots(IPA1Toolbox.species);
		Map<Gene,Set<RecS>> genesWithRef = getAllGenesWithRefMotif(genes,IPA1Toolbox.reference);
		
		System.out.println("Number of osa genes in Monocots: "+ genes.size());
		System.out.println("Number of osa genes containing reference motif: " + genesWithRef.size());
				
		
		IPA1Toolbox.printTableS(genesWithRef, IPA1Toolbox.getFilenamePWMMatches());

		
		String v1="TGGGCC";
		Map<Gene,Set<RecS>> genesWithV1 = getAllGenesWithRefMotif(genes,v1);
		
		System.out.println("Number of osa genes in Monocots: "+ genes.size());
		System.out.println("Number of osa genes containing "+v1+" : " + genesWithV1.size());
		

		String v2="TGGGCT";
		Map<Gene,Set<RecS>> genesWithV2 = getAllGenesWithRefMotif(genes,v2);
		
		System.out.println("Number of osa genes in Monocots: "+ genes.size());
		System.out.println("Number of osa genes containing "+v2+" : " + genesWithV2.size());
		
		
	}
}




