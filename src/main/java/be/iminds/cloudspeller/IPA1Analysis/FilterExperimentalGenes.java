package be.iminds.cloudspeller.IPA1Analysis;

import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.Sequence;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.KN1Analysis.RecF;
import be.iminds.cloudspeller.KN1Analysis.RecS;

public class FilterExperimentalGenes {


	
	
	
	public static void main (String [] args) throws IOException{
		
		filterExperimentalIPA1Genes();
		
	}
	
	
	public static void filterExperimentalIPA1Genes() throws IOException{
		String species = "OS";
		
		//different be.iminds.cloudspeller.input files: IPA1_BR -> bound + regulated
		//						 IPA1_YP
		//						 IPA1_SA
		//						 IPA1_Union
		//String data = IPA1_Inters
		
		String data = "IPA1_BR";
		String inputFile = IPA1Toolbox.getIPA1Directory()+data;
		@SuppressWarnings("unused")
		Set<Gene> IPAGenesBR = getGenesFromFile(inputFile,species);

		data = "IPA1_YP";
		inputFile = IPA1Toolbox.getIPA1Directory()+data;
		Set<Gene> IPAGenesYP = getGenesFromFile(inputFile,species);
		
		data = "IPA1_SA";
		inputFile = IPA1Toolbox.getIPA1Directory()+data;
		Set<Gene> IPAGenesSA = getGenesFromFile(inputFile,species);
		
		Set<Gene> IPAGenesUnion = new HashSet<Gene>(); 
		IPAGenesUnion.addAll(IPAGenesSA); IPAGenesUnion.addAll(IPAGenesYP);
		
		Set<Gene> IPAGenesIntersection =IPA1Toolbox.calculateIntersection(IPAGenesYP,IPAGenesSA); 
		
		Set<Gene> IPAGenes = IPAGenesIntersection;
		
				
		IPA1Toolbox.setDatasetFromFile();
		
		Map<Gene,Set<String>> IPAFamilyMap = getFamilyOfGenes(IPAGenes,species);
		
		System.out.println("Stats:");
		System.out.println("Number of IPA1 genes in paper: "+ IPAGenes.size());
		System.out.println("Number of IPA1 genes in Monocot dataset: "+ IPAFamilyMap.size());
		
		Map<Gene, Set<RecS>> tableSmall = getRecordsWithReferencesInPromoter(IPAFamilyMap);
		System.out.println("TABLESMALL #query genes: "+tableSmall.size());
		
		IPA1Toolbox.printTableS(tableSmall, IPA1Toolbox.getFilenameExperimentalMatches());

		
		Map<Gene, Set<RecF>> tableFull  = getRecordsWithReferenceInPromoterAndOrthologs(IPAFamilyMap,species);
		System.out.println("TABLEFULL #query genes: "+tableFull.size());

		IPA1Toolbox.printTableF(tableFull, IPA1Toolbox.getFilenameExperimentalMatchesWithOrthologs());

		
	}


	private static Map<Gene, Set<RecF>> getRecordsWithReferenceInPromoterAndOrthologs(
			Map<Gene, Set<String>> iPAFamilyMap, String species) {

		Map<Gene, Set<RecF>> table = new HashMap<Gene, Set<RecF>>();
		
		for (Map.Entry<Gene, Set<String>> e : iPAFamilyMap.entrySet()){
			
			for (String gfid : e.getValue()){
			
				GeneFamily gf = IPA1Toolbox.getFamilies().get(gfid);
				ArrayList<Gene> genes = gf.getGenes();
				ArrayList<Sequence> sequences = gf.getSequences();
				
				Sequence sequence = null;
				for (int i=0; i<genes.size(); i++){
					if (genes.get(i).equals(e.getKey())){ 
						sequence = sequences.get(i);
						break;
					}
				}
					
				List<Suffix> matches = IPA1Toolbox.giveMatchesOfPatternInSequence(IPA1Toolbox.reference,sequence);
							
				if (matches!=null){ //if query gene has ref motif in promoter check orthologs
				
					Set<RecF> records = new HashSet<RecF>();
					
					Set<String> orthos = new HashSet<String>();
					for (int i=0; i<genes.size(); i++){
						Sequence currentSeq = sequences.get(i);
						
						List<Suffix> tempMatches = IPA1Toolbox.giveMatchesOfPatternInSequence(IPA1Toolbox.reference,currentSeq);
						
						if (tempMatches!=null){
							
							if (!genes.get(i).getOrganism().equals(species)){ //check if match in ortholog => records will be added
								orthos.add(genes.get(i).getOrganism());
							}
							
							for (Suffix s : tempMatches){
								String bs = IPA1Toolbox.getActualBS(s,currentSeq);
								records.add(new RecF(gfid,genes.get(i).getID(),bs,-1));
							}
						}
					}
				
					//debug orthos.size==1
					if (orthos.size()>0 /*==3*/ ){
						table.put(e.getKey(),records);
					}
				
				} else {
				//No matches with ref motif, do nothing!
				}
			}
		}
		return table;
	}


	private static Map<Gene, Set<RecS>> getRecordsWithReferencesInPromoter(
			Map<Gene, Set<String>> iPAFamilyMap) {
	
		Map<Gene, Set<RecS>> table = new HashMap<Gene, Set<RecS>>();
		
		for (Map.Entry<Gene, Set<String>> e : iPAFamilyMap.entrySet()){

			String gfid = e.getValue().iterator().next();
			GeneFamily gf = IPA1Toolbox.getFamilies().get(gfid);
			
			ArrayList<Gene> genes = gf.getGenes();
			ArrayList<Sequence> sequences = gf.getSequences();
								
			Set<RecS> records = new HashSet<RecS>(); 
							
			Sequence sequence = null;
			for (int i=0; i<genes.size(); i++){
				if (genes.get(i).equals(e.getKey())){ //is query gene?
					sequence = sequences.get(i);
					break;
				}
			}
			
			List<Suffix> matches = IPA1Toolbox.giveMatchesOfPatternInSequence(IPA1Toolbox.reference,sequence);
			
			if (matches!=null){
			
				for (Suffix s : matches){
					String bs = IPA1Toolbox.getActualBS(s, sequence);
					records.add(new RecS(bs,-1));
				}
				
				table.put(e.getKey(),records);
			} else {
				//No matches with ref motif, do nothing!
			}
		}
		return table;
	}

		
	


	private static Map<Gene, Set<String>> getFamilyOfGenes(Set<Gene> genes, String species) {
	
		Map<String,GeneFamily> families = IPA1Toolbox.getFamilies();
		Map<Gene, Set<String>> geneToOrthoMap = new HashMap<Gene, Set<String>>();
		
		for (Map.Entry<String,GeneFamily> gfp : families.entrySet()){
			for (Gene g : gfp.getValue().getGenes()){
				if (g.getOrganism().equals(species) && genes.contains(g)){ 
					
					Set<String> familiesContainingThisGene = geneToOrthoMap.get(g);
					if (familiesContainingThisGene==null){
						familiesContainingThisGene = new HashSet<String>();
						geneToOrthoMap.put(g,familiesContainingThisGene);
					}
					
					familiesContainingThisGene.add(gfp.getValue().getFamilyName());
				}
			}
		}
		
		return geneToOrthoMap;
	
	}
		
		
	


	/**
	 * @return The list of osa genes regulated by the IPA1 TF
	 */
	public static Set<Gene> getGenesFromFile(String inputFile, String species) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		
		Set<Gene> geneIDs = new HashSet<Gene>();
		String line;
		while ((line=in.readLine())!=null){
			geneIDs.add(new Gene(line.trim(),species)); //gene id
		}
		in.close();
		return geneIDs;
	}
	


}
