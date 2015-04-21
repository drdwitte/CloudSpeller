package be.iminds.cloudspeller.KN1Analysis;

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


public class PWMScanDenovo {
	
	Map<Gene, Set<RecS>> getAllMaizeGenesWithRefMotif(
			Map<Gene, Sequence> maizeGenes, String referenceMotif) {
		int c=0;
		Map<Gene, Set<RecS>> recordMap = new HashMap<Gene, Set<RecS>>();
		
		for (Map.Entry<Gene,Sequence> e : maizeGenes.entrySet()){
			c++;
			if (c%100==0){
				System.out.println("Maize promoters scanned: "+c);
			}
			Sequence seq = e.getValue();
			List<Suffix> matches = KN1Toolbox.giveMatchesOfPatternInSequence(referenceMotif, seq);
						
			if (matches!=null){
				Set<RecS> records = new HashSet<RecS>();
				for (Suffix m : matches){
					String bs = KN1Toolbox.getActualBS(m, seq);
					records.add(new RecS(bs,-1));
				}
				
				recordMap.put(e.getKey(),records);
			}
		
		}
		return recordMap;
		
	}

	Map<Gene, Sequence> getAllMaizeGeneDataFromMonocots() {
		
		Map<Gene, Sequence> geneSeqs = new HashMap<Gene, Sequence>();
		
		for (Map.Entry<String,GeneFamily> e : KN1Toolbox.getFamilies().entrySet()){
			
			ArrayList<Gene> genes = e.getValue().getGenes();
			ArrayList<Sequence> seqs = e.getValue().getSequences();
			
			for (int i=0; i<genes.size(); i++){
				Gene g = genes.get(i);
				if (g.getOrganism().equals("ZM")){ //is maize gene?
					
					geneSeqs.put(genes.get(i),seqs.get(i));
				}
			}
		} 
		return geneSeqs;
	}
	
	public static void main (String [] args) throws IOException{
		
		String filenameDataset = KN1Toolbox.getDatasetPath();
		PWMScanDenovo scanner = new PWMScanDenovo();
		KN1Toolbox.setDatasetFromFile(filenameDataset);
						
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		
		Map<Gene,Sequence> maizeGenes = scanner.getAllMaizeGeneDataFromMonocots();
		Map<Gene,Set<RecS>> maizeGenesWithRef = scanner.getAllMaizeGenesWithRefMotif(maizeGenes,KN1Toolbox.referenceMotif);
		
		System.out.println("Number of maize genes in Monocots: "+ maizeGenes.size());
		System.out.println("Number of maize genes containing reference motif: " + maizeGenesWithRef.size());
				
		KN1Toolbox.generateWeightForAllMatches(maizeGenesWithRef,null,pwm);
			
		//Filter on PWM score
		Map<Gene, Set<RecS>> maizeGenesFiltered;
		
		int t = KN1Toolbox.getPWMScoreThreshold();
		String filename = KN1Toolbox.generateOutputFilenamePWMScan(t);
		
		
		maizeGenesFiltered = KN1Toolbox.selectRecSAbovePWMThreshold(maizeGenesWithRef,t);
		System.out.println("#monocot maize genes with filter = "+t +" : "+maizeGenesFiltered.size());
		KN1Toolbox.printTableS(maizeGenesFiltered, filename);
		
	}
}
