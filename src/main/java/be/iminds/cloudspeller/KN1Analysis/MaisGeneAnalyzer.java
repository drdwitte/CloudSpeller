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


public class MaisGeneAnalyzer {

	public static final double lowerBoundScore = -1000000;
	private static Map<String,String> degMatchMap;
    static {
		degMatchMap = new HashMap<String, String>();
		degMatchMap.put("AA","A");
		degMatchMap.put("AC","M");
		degMatchMap.put("AG","R");
		degMatchMap.put("AT","W");
		degMatchMap.put("CA","M");
		degMatchMap.put("CC","C");
		degMatchMap.put("CG","S");
		degMatchMap.put("CT","Y");
		degMatchMap.put("GA","R");
		degMatchMap.put("GC","S");
		degMatchMap.put("GG","G");
		degMatchMap.put("GT","K");
		degMatchMap.put("TA","W");
		degMatchMap.put("TC","Y");
		degMatchMap.put("TG","K");
		degMatchMap.put("TT","T");
    }
	

	public static void main (String [] args) throws IOException{
		
		String bolducFile = args[0];
		String filenameDataset = args[1];
		
		PWMScanDenovo scanner = new PWMScanDenovo();
		BolducAnalyzer bolducAnalyzer = new BolducAnalyzer();
		KN1Toolbox.setDatasetFromFile(filenameDataset);
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		Set<Gene> bolducGenes = bolducAnalyzer.getBolducGenesFromFile(bolducFile);
		
		Map<Gene,Sequence> maizeGenes = scanner.getAllMaizeGeneDataFromMonocots();
		Map<Gene,Set<RecS>> maizeGenesWithRef = scanner.getAllMaizeGenesWithRefMotif(maizeGenes,KN1Toolbox.referenceMotif);
		
		Map<Gene,Set<String>> maizeGeneFamilyMap = bolducAnalyzer.getFamilyOfGenes(maizeGenes.keySet());

		//add pwmscores for binding sites in maizeGene
		KN1Toolbox.generateWeightForAllMatches(maizeGenesWithRef,null,pwm); //generate pwm scores
		
		for (Map.Entry<Gene,Set<RecS>> e : maizeGenesWithRef.entrySet()){
			
			Gene maizeGene = e.getKey();
			Set<RecS> records = e.getValue();
			
			StringBuilder bsInMaize = new StringBuilder();
			Set<String> bindingSitesInMaizeQueryGene = new HashSet<String>();
					
			double maxMaizeScore = lowerBoundScore;
			
			for (RecS record : records){
							
				bsInMaize.append(record.getBindingSite());
				bsInMaize.append("\t");
				double score = record.getPWMScore();
				bsInMaize.append(score);
				bsInMaize.append("\n");
				bindingSitesInMaizeQueryGene.add(record.getBindingSite());
			

				if (score > maxMaizeScore){
					maxMaizeScore = score;
				}
			}
			
			
			
			
			boolean isBolduc = bolducGenes.contains(maizeGene);
			char bolduc = isBolduc?'y':'n';
			System.out.println("*");
			System.out.println(maizeGene.getID().toString()+"\t"+bolduc);
			System.out.println("BS in maize Gene:");
			System.out.print(bsInMaize);
			System.out.println("MaxMaizeScore: "+maxMaizeScore);
			
			Set<String> familyOccs = maizeGeneFamilyMap.get(maizeGene);
			System.out.println("GeneFamilies: "+familyOccs.size());
			
			double maxOrthoScore = lowerBoundScore;
			Set<String> allBSitesInOrthologs = new HashSet<String>();
			
			
			for (String f : familyOccs){
				
				System.out.println(f); //print gene families in which gene occurs
				
				GeneFamily gf = KN1Toolbox.getFamilies().get(f);
				ArrayList<Gene> genes = gf.getGenes();
				ArrayList<Sequence> sequences = gf.getSequences();
				
				
				
				//all binding sites (pwmscore currently -1)
				Set<RecF> newRecordsForOrtho = new HashSet<RecF>();
				
				
				for (int i=0; i<genes.size(); i++){
					
					if (genes.get(i).getID().startsWith("ZM")){
						continue; //skip query maizegene and paralogs
					}
					
					Sequence sequence = sequences.get(i);
					List<Suffix> matches = KN1Toolbox.giveMatchesOfPatternInSequence(KN1Toolbox.referenceMotif,sequence);

					if (matches==null){
						continue;
					}
										
					for (Suffix s : matches){
						String bs = KN1Toolbox.getActualBS(s,sequence);
						newRecordsForOrtho.add(new RecF(f,genes.get(i).getID(),bs,-1));
						allBSitesInOrthologs.add(bs);
					}
				}
				
				//map required to comply with signature of pwmscore calculation
				Map<Gene,Set<RecF>> orthoRecords = new HashMap<Gene,Set<RecF>>();
				orthoRecords.put(maizeGene,newRecordsForOrtho);

				//generate pwm score for ortholog binding sites
				KN1Toolbox.generateWeightForAllMatches(null,orthoRecords,pwm);
								
				//print all bindingsites with pwmscore
				for (RecF r : orthoRecords.get(maizeGene)){	
					double score = r.getPWMScore();
					System.out.println(r.getGene()+"\t"+r.getBindingSite()+"\t"+score);
					if (score > maxOrthoScore){
						maxOrthoScore = score;
					}
				}
			}
			
			System.out.println("HighestOrthologScore: "+ maxOrthoScore);
			
			Set<String> matchingMotifs = getLDMotifsMatchingInMaizeAndOrtholog(bindingSitesInMaizeQueryGene,allBSitesInOrthologs);
			
			System.out.println("Least Degenerate Matching motifs (multispecies match):");
			for (String matchingMotif : matchingMotifs){
				System.out.print(matchingMotif+" ");
			}
			System.out.print("\n");
			
			bsInMaize.setLength(0);
			System.out.println("");
			
			
		}

		
	}

	private static Set<String> getLDMotifsMatchingInMaizeAndOrtholog(
			Set<String> bindingSitesInMaizeQueryGene,
			Set<String> allBSitesInOrthologs) {
		
		Set<String> allMatches = new HashSet<String>();
		for (String s1 : bindingSitesInMaizeQueryGene){
			for (String s2 : allBSitesInOrthologs){
				allMatches.add(LDMotifIntersection(s1,s2));
			}
		}
		
		return allMatches;
	}

	public static String LDMotifIntersection(String s1,
			String s2) {
	
		String regex = KN1Toolbox.regex;
	
		List<StringBuilder> oldA = new ArrayList<StringBuilder>();
		List<StringBuilder> newA = new ArrayList<StringBuilder>();
		
		oldA.add(new StringBuilder(""));
		
		for (int i=0; i<regex.length(); i++){
			char c = regex.charAt(i);
			
			if (c == '.'){
				
				String charInBothStrings = ""+s1.charAt(i)+s2.charAt(i);
				String alph = leastDegMatches(charInBothStrings);			
				
				for (StringBuilder b : oldA){
					
					for (int j=0; j<alph.length(); j++){
						
						StringBuilder copy = new StringBuilder(b);
						copy.append(alph.charAt(j));
						newA.add(copy);
					}
				}
				oldA.clear();
				oldA.addAll(newA);
				
				newA.clear();
				//oldA = newA;
				
			} else {
				for (StringBuilder b : oldA){
					b.append(c);
				}
			}
		}
		
		
		/*Set<String> intersection = new HashSet<String>();
		
		for (StringBuilder sb : oldA){
			intersection.add(sb.toString());
		}*/
		
		return oldA.get(0).toString();
	}

	public static String leastDegMatches(String charInBothStrings) {
		return degMatchMap.get(charInBothStrings);
	}

	
	
	
	

	
	
	
	
	
	
}
