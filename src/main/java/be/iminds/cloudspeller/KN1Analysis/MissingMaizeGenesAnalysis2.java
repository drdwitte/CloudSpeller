package be.iminds.cloudspeller.KN1Analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.input.Gene;

public class MissingMaizeGenesAnalysis2 {

	
	public static void main(String [] args) throws IOException{
		
		Gene [] maizeGenesMissing = {
			new Gene("ZM02G13060","ZM"),
			new Gene("ZM09G06540","ZM"), 
			new Gene("ZM02G38700","ZM"),
			new Gene("ZM05G37850","ZM"), 
			new Gene("ZM06G24040","ZM")
		};
		
		Set<Gene> missingGenes = new HashSet<Gene>();
		for (Gene g : maizeGenesMissing){
			missingGenes.add(g);
		}
		
		String dataset = KN1Toolbox.getDatasetPath();
		
		BolducAnalyzer bolducAnalyzer = new BolducAnalyzer();
		KN1Toolbox.setDatasetFromFile(dataset);
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		Map<Gene,Set<String>> bolducFamilyMap = bolducAnalyzer.getFamilyOfGenes(missingGenes);

		System.out.println(bolducFamilyMap.size());
		
		Map<Gene, Set<RecF>> tableFull  = bolducAnalyzer.getBolducRecordsWithReferenceInPromoterAndOrthologs(bolducFamilyMap,KN1Toolbox.referenceMotif);
		System.out.println("TABLEFULL #query genes: "+tableFull.size());
		
		System.out.println("Adding scores...");
		KN1Toolbox.generateWeightForAllMatches(null,tableFull,pwm);
		System.out.println("done");
		
		
		DegMatchesProcessor processor = new DegMatchesProcessor(null);
		processor.initializeBSMap(KN1Toolbox.regex, pwm);
		
		
		for (Map.Entry<Gene, Set<RecF>> e : tableFull.entrySet()){
			System.out.println(e.getKey());
			
			//new motif moet met ortholog dus geen maize!!
			Set<RecF> queryRec = new HashSet<RecF>();
			Set<RecF> orthologs = new HashSet<RecF>();
			
			for (RecF r : e.getValue()){
				
				if (r.getGene().equals(e.getKey().getID())){
					if (r.getPWMScore()>=0.0){
						queryRec.add(r);
					}
				} else if (r.getGene().startsWith("ZM")){
					//ignore
				} else {
					
					if (r.getPWMScore()>=0.0){
						orthologs.add(r);
				
					}
				}
			}
		
			for (RecF r1 : queryRec)
				for (RecF r2 : orthologs){
				String newMot = LDMotifIntersection(r1.getBindingSite(),r2.getBindingSite());
				
				System.out.println(r2+"\t->"+newMot);
				System.out.println("Stats: pwm= "+processor.calculatescoreWorstBS(newMot));
				System.out.println("Number of deg positions"+numDegPos(newMot));
				
			}
			
			
			System.out.println("");
		}
		
		
		//get all orthologous records
		
		//add pwm scores
		
		//print stats about why maybe missed: alphabet + number of deg positions
		
	}
	
	private static int numDegPos(String w) {
		
		int n=0;
		
		for (int i=0; i<w.length(); i++){
			
			switch (w.charAt(i)){
			case 'A': break;
			case 'C': break;
			case 'G': break;
			case 'T': break;
			default:
				n++;
			}
			
		}
		return n;
		
		
	}

	/*private static Set<String> getLDMotifsMatchingInMaizeAndOrtholog(
			Set<String> bindingSitesInMaizeQueryGene,
			Set<String> allBSitesInOrthologs) {
		
		Set<String> allMatches = new HashSet<String>();
		for (String s1 : bindingSitesInMaizeQueryGene){
			for (String s2 : allBSitesInOrthologs){
				allMatches.add(LDMotifIntersection(s1,s2));
			}
		}
		
		return allMatches;
	}*/
	
	public static String LDMotifIntersection(String s1,
			String s2) {
	
		String regex = KN1Toolbox.regex;

		String degMotif="";
		
		for (int i=0; i<regex.length(); i++){
			char c = regex.charAt(i);
			
			if (c == '.'){
				
				String charInBothStrings = ""+s1.charAt(i)+s2.charAt(i);
				String alph = leastDegMatches(charInBothStrings);			
				
				degMotif+=alph;
				
				
			} else {
				degMotif+=c;
			}
		}
		return degMotif;
	}
	
	public static String leastDegMatches(String charInBothStrings) {
		return degMatchMap.get(charInBothStrings);
	}
	
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
}
