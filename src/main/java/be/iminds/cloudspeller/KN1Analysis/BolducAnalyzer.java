package be.iminds.cloudspeller.KN1Analysis;

import be.iminds.cloudspeller.indexing.FullMotifMatchWithPos;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.indexing.SuffixInterpreter;
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
import java.util.Scanner;
import java.util.Set;



/**
 * @author ddewitte
 */

public class BolducAnalyzer {

	public BolducAnalyzer(){}
	
	/**
	 * @return The list of genes regulated by the KN1 TF according to Bolduc et al
	 */
	public Set<Gene> getBolducGenesFromFile(String filename) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(filename));
		
		Set<Gene> geneIDs = new HashSet<Gene>();
		String line;
		while ((line=in.readLine())!=null){
			Scanner scanner = new Scanner(line);
			scanner.next(); //genome id
			geneIDs.add(new Gene(scanner.next(),"ZM")); //gene id
		}
		in.close();
		return geneIDs;
	}
	
	/**
	 * @return Mapping between Mais Query gene and orthologous gene families in Monocot data
	 */
	public Map<Gene, Set<String>> getFamilyOfGenes(Set<Gene> bolducGenes){
	
		Map<String,GeneFamily> families = KN1Toolbox.getFamilies();
		Map<Gene, Set<String>> bolducToOrthoMap = new HashMap<Gene, Set<String>>();
		
		for (Map.Entry<String,GeneFamily> gfp : families.entrySet()){
			for (Gene g : gfp.getValue().getGenes()){
				if (g.getOrganism().equals("ZM") && bolducGenes.contains(g)){ //maizeGene
					
					Set<String> familiesContainingThisGene = bolducToOrthoMap.get(g);
					if (familiesContainingThisGene==null){
						familiesContainingThisGene = new HashSet<String>();
						bolducToOrthoMap.put(g,familiesContainingThisGene);
					}
					
					familiesContainingThisGene.add(gfp.getValue().getFamilyName());
				}
			}
		}
		
		return bolducToOrthoMap;
	
	}
	
	/**
	 * @return Table with structure: Mais Query gene - Binding site - PWM score for BS
	 * One line per KN1 binding site (note multiple bs per sequence possible!)
	 */
	public Map<Gene, Set<RecS>> getBolducRecordsWithReferenceInPromoter(Map<Gene,Set<String>> bolducFamilyMapping, String reference){
		int c=0;		
		Map<Gene, Set<RecS>> table = new HashMap<Gene, Set<RecS>>();
		
		for (Map.Entry<Gene, Set<String>> e : bolducFamilyMapping.entrySet()){
			
			c++;
			if (c%100==0){
				System.out.println("Number of query genes processed: "+c);
			}
			
			String gfid = e.getValue().iterator().next();
			GeneFamily gf = KN1Toolbox.getFamilies().get(gfid);
			
			ArrayList<Gene> genes = gf.getGenes();
			ArrayList<Sequence> sequences = gf.getSequences();
								
			Set<RecS> records = new HashSet<RecS>(); 
							
			Sequence maizeSequence = null;
			for (int i=0; i<genes.size(); i++){
				if (genes.get(i).equals(e.getKey())){ //is bolduc query gene?
					maizeSequence = sequences.get(i);
					break;
				}
			}
			
			List<Suffix> matches = KN1Toolbox.giveMatchesOfPatternInSequence(reference,maizeSequence);
			
			if (matches!=null){
			
				for (Suffix s : matches){
					String bs = KN1Toolbox.getActualBS(s, maizeSequence);
					records.add(new RecS(bs,-1));
				}
				
				table.put(e.getKey(),records);
			} else {
				//No matches with ref motif, do nothing!
			}
		}
		return table;
	}

	/**
	 * @return Table with structure: Mais Query gene - GeneFamily - Gene - Binding site - PWM score for BS
	 * One line per KN1 binding site in any of the bolduc orthologs
	 */
	public Map<Gene, Set<RecF>> getBolducRecordsWithReferenceInPromoterAndOrthologs(Map<Gene,Set<String>> bolducFamilyMapping, String reference){

		int c=0;
		Map<Gene, Set<RecF>> table = new HashMap<Gene, Set<RecF>>();
		
		for (Map.Entry<Gene, Set<String>> e : bolducFamilyMapping.entrySet()){
			c++;
			if (c%100==0){
				System.out.println("Number of query genes processed: "+c);
			}
			
			for (String gfid : e.getValue()){
			
				GeneFamily gf = KN1Toolbox.getFamilies().get(gfid);
				ArrayList<Gene> genes = gf.getGenes();
				ArrayList<Sequence> sequences = gf.getSequences();
				
				//find the bolduc query gene
				Sequence maizeSequence = null;
				for (int i=0; i<genes.size(); i++){
					if (genes.get(i).equals(e.getKey())){ 
						maizeSequence = sequences.get(i);
						break;
					}
				}
					
				List<Suffix> matches = KN1Toolbox.giveMatchesOfPatternInSequence(reference,maizeSequence);
							
				if (matches!=null){ //if query gene has ref motif in promoter check orthologs
				
					Set<RecF> records = new HashSet<RecF>();
					boolean matchInOrtholog = false;
					for (int i=0; i<genes.size(); i++){
						Sequence currentSeq = sequences.get(i);
						
						List<Suffix> tempMatches = KN1Toolbox.giveMatchesOfPatternInSequence(reference,currentSeq);
						
						if (tempMatches!=null){
							
							if (!genes.get(i).getOrganism().equals("ZM")){ //check if match in ortholog => records will be added
								matchInOrtholog=true;
							}
							
							for (Suffix s : tempMatches){
								String bs = KN1Toolbox.getActualBS(s,currentSeq);
								records.add(new RecF(gfid,genes.get(i).getID(),bs,-1));
							}
						}
					}
				
				
					if (matchInOrtholog){
						table.put(e.getKey(),records);
					}
				
				} else {
				//No matches with ref motif, do nothing!
				}
			}
		}
		return table;
	}
	
	public Map<Gene, Set<RecF>> selectBRecFAbovePWMThreshold(Map<Gene, Set<RecF>> tableFull, int pwmThreshold){
		
		Map<Gene, Set<RecF>> filteredTable = new HashMap<Gene, Set<RecF>>();
		
		for (Map.Entry<Gene,Set<RecF>> e : tableFull.entrySet()){
			
			
			Map<Gene,Set<RecF>> filteredRecords = new HashMap<Gene,Set<RecF>>();
			
			for (RecF unfilteredRec : e.getValue()){ //eliminate all records which do not satisfy pwm threshold
				
				if (unfilteredRec.getPWMScore() >=pwmThreshold){
					
					Set<RecF> tempFilteredRecs = filteredRecords.get(e.getKey());
					
					if (tempFilteredRecs==null){
						tempFilteredRecs = new HashSet<RecF>();
						filteredRecords.put(e.getKey(),tempFilteredRecs);
					}
					tempFilteredRecs.add(unfilteredRec);
				}
			}
			
			
			
			//check if match with maize query gene left and with at least one ortholog for every genefamily 
			// in which the query gene occurs
			for (Map.Entry<Gene,Set<RecF>> p : filteredRecords.entrySet()){
				
				boolean matchWithQueryG = false;
				boolean matchWithOrtholog = false;
				
				for (RecF rec : p.getValue()){
					if (!rec.getGene().startsWith("ZM")){
						matchWithOrtholog = true;
					} else if (rec.getGene().equals(p.getKey().getID())){
						matchWithQueryG = true;
					}
				}
				
				
				if (matchWithOrtholog && matchWithQueryG){
					
					
					Set<RecF> filteredRecsAlreadyStored = filteredTable.get(e.getKey());
					
					if (filteredRecsAlreadyStored == null){
						filteredRecsAlreadyStored = new HashSet<RecF>();
						filteredTable.put(e.getKey(),filteredRecsAlreadyStored);
					}
					filteredRecsAlreadyStored.addAll(p.getValue());
				}
			}
		}
		return filteredTable;
		
	}
	
	
	public static void filterBolducGenes() throws IOException{
		
		String bolducFile = KN1Toolbox.getKN1Directory()+"KN1TargetsAllFromBolduc.txt";
		String dataset = KN1Toolbox.getDatasetPath();
				
		BolducAnalyzer bolducAnalyzer = new BolducAnalyzer();
		KN1Toolbox.setDatasetFromFile(dataset);
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		Set<Gene> bolducGenes = bolducAnalyzer.getBolducGenesFromFile(bolducFile);
		Map<Gene,Set<String>> bolducFamilyMap = bolducAnalyzer.getFamilyOfGenes(bolducGenes);
		
		System.out.println("Stats:");
		System.out.println("Number of bolduc genes in KN1 paper: "+ bolducGenes.size());
		System.out.println("Number of bolduc genes in Monocot dataset: "+ bolducFamilyMap.size());
		
		Map<Gene, Set<RecS>> tableSmall = bolducAnalyzer.getBolducRecordsWithReferenceInPromoter(bolducFamilyMap,KN1Toolbox.referenceMotif);
		System.out.println("TABLESMALL #query genes: "+tableSmall.size());
		
		Map<Gene, Set<RecF>> tableFull  = bolducAnalyzer.getBolducRecordsWithReferenceInPromoterAndOrthologs(bolducFamilyMap,KN1Toolbox.referenceMotif);
		System.out.println("TABLEFULL #query genes: "+tableFull.size());
		
		System.out.println("Adding scores...");
		KN1Toolbox.generateWeightForAllMatches(tableSmall,tableFull,pwm);
		System.out.println("done");
		
		//Filter on PWM score
		Map<Gene, Set<RecS>> tableSmallFiltered;
		
		int t = KN1Toolbox.getPWMScoreThreshold();
		String filename = KN1Toolbox.generateOutputFilenameBolduc(t);

		tableSmallFiltered = KN1Toolbox.selectRecSAbovePWMThreshold(tableSmall,t);
		System.out.println("#Bolduc query genes with filter = "+t +" : "+tableSmallFiltered.size());
		KN1Toolbox.printTableS(tableSmallFiltered, filename);
		
		Map<Gene, Set<RecF>> tableFullFiltered;
		filename = KN1Toolbox.generateOutputFilenameBolducWithOrthologs(t);
	
			
		System.out.println("Threshold = "+t);
	
		tableFullFiltered = bolducAnalyzer.selectBRecFAbovePWMThreshold(tableFull,t);
		System.out.println("#BolducO query genes with filter = "+t +" : "+tableFullFiltered.size());
		KN1Toolbox.printTableF(tableFullFiltered, filename);
		
		
		//NOTE genen uit records table halen is gewoon .keySet() gebruiken!

	}

	public static void findBolducTargets() throws IOException{
		
		String bolducFile = KN1Toolbox.getKN1Directory()+"KN1TargetsAllFromBolduc.txt";
		String dataset = KN1Toolbox.getDatasetPath();
				
		BolducAnalyzer bolducAnalyzer = new BolducAnalyzer();
		KN1Toolbox.setDatasetFromFile(dataset);
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		Set<Gene> bolducGenes = bolducAnalyzer.getBolducGenesFromFile(bolducFile);
		Map<Gene,Set<String>> bolducFamilyMap = bolducAnalyzer.getFamilyOfGenes(bolducGenes);
		
		System.out.println("Stats:");
		System.out.println("Number of bolduc genes in KN1 paper: "+ bolducGenes.size());
		System.out.println("Number of bolduc genes in Monocot dataset: "+ bolducFamilyMap.size());
		
		Map<Gene, Set<RecS>> tableSmall = bolducAnalyzer.getBolducRecordsWithReferenceInPromoter(bolducFamilyMap,KN1Toolbox.referenceMotif);
		System.out.println("TABLESMALL #query genes: "+tableSmall.size());
		
		System.out.println("Adding scores...");
		KN1Toolbox.generateWeightForAllMatches(tableSmall,null,pwm);
		System.out.println("done");
		
		//Filter on PWM score
		Map<Gene, Set<RecS>> tableSmallFiltered;
		
		int t = KN1Toolbox.getPWMScoreThreshold();

		tableSmallFiltered = KN1Toolbox.selectRecSAbovePWMThreshold(tableSmall,t);
		System.out.println("#Bolduc query genes with filter = "+t +" : "+tableSmallFiltered.size());
		
			//GET TARGETS
		
		for (Map.Entry<Gene, Set<RecS>> e : tableSmallFiltered.entrySet()){
			
			for (RecS record : e.getValue()){
				String bs = record.getBindingSite();
				Gene g = e.getKey();
				
				String targetFamilyID = bolducFamilyMap.get(e.getKey()).iterator().next();
				GeneFamily gf = KN1Toolbox.getFamilies().get(targetFamilyID);
				
				
				Sequence maizeSequence = gf.getSequence(g);
				List<Suffix> matches = KN1Toolbox.giveMatchesOfPatternInSequence(bs,maizeSequence);
				
				for (int i=0; i<matches.size(); i++){
				
					FullMotifMatchWithPos outputLine =  SuffixInterpreter.translateSuffixWithPosSingleSeq
						(bs,matches.get(i), targetFamilyID, g.getID(), maizeSequence.length()); 
					
					System.out.println(outputLine);
				}
			}
			
			
		}
	}
	
	
	public static void main (String [] args) throws IOException {
		//filterBolducGenes();
		
		findBolducTargets();
				
	}
}
