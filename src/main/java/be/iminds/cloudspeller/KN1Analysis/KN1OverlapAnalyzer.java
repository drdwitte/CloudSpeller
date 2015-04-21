package be.iminds.cloudspeller.KN1Analysis;

import be.iminds.cloudspeller.input.Gene;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;



public class KN1OverlapAnalyzer {

	
	private void singleMotifOverlap(String exp) throws IOException {
		
		int t = KN1Toolbox.getPWMScoreThreshold();
		Map<Gene, Set<RecS>> discTableA;
		discTableA = KN1Toolbox.readTableSFromFile(KN1Toolbox.generateOutputFilenameDeNovoWorstBS(exp,t));
		
		//read DeNovo records and put them in a sortedMap (scorePair)
		String outputRewWorstScore = KN1Toolbox.getDeNovoRecords(exp);
		
		Set<DeNovoRecord> kn1Variants = DegMatchesProcessor.selectDeNovoRecords(t,outputRewWorstScore);
		SortedMap<ScorePair,DeNovoRecord> motifPairsMap = generateMotifMapFromDeNovoRecords(kn1Variants);// new TreeMap<ScorePair,DeNovoRecord>();
		
		//OVERLAP WITH BOLDUC
		
		//read bolduc table for pwmthreshold 0
		System.out.println("Reading: "+KN1Toolbox.generateOutputFilenameBolduc(t));
		Map<Gene, Set<RecS>> bolducTable0 = KN1Toolbox.readTableSFromFile(KN1Toolbox.generateOutputFilenameBolduc(t));

		
		System.out.println("Number of maizeGenes in Bolduc with threshold "+t+" : "+bolducTable0.keySet().size());
		
		Set<Gene> Bolduc0 = new HashSet<Gene>(bolducTable0.keySet());
		//Set<Gene> BolducA = new HashSet<Gene>(bolducTableA.keySet());
				
		
		System.out.println("variant \t blsMin \t #Fam \t Conf \t score \t #maize targets \t  ovBold \t Pwmthreshold: "+t);
		performOverlapAnalysis(motifPairsMap,Bolduc0,discTableA);
		
		System.out.println(""); System.out.println(""); System.out.println("");
		
		//OVERLAP WITH BOLDUCO
	
		
		//read bolduc with orthologs table for pwmthreshold 0 and -1000 (=All)
		Map<Gene, Set<RecF>> bolducOTable0 = KN1Toolbox.readTableFFromFile(KN1Toolbox.generateOutputFilenameBolducWithOrthologs(t));
		//Map<Gene, Set<RecF>> bolducOTableA = KN1Toolbox.readTableFFromFile(KN1Toolbox.generateOutputFilenamesBolducWithOrthologs().get(0));
		
		System.out.println("Number of maizeGenes in BoldO with threshold "+t+" : "+bolducOTable0.keySet().size());

		Set<Gene> BolducO0 = new HashSet<Gene>(bolducOTable0.keySet());
		//Set<Gene> BolducOA = new HashSet<Gene>(bolducOTableA.keySet());

		System.out.println("variant \t blsMin \t #Fam \t Conf \t score \t #maize targets \t ovBold(0)");
		performOverlapAnalysis(motifPairsMap,BolducO0,discTableA);
		
		
	}
	
	private void performOverlapAnalysis(SortedMap<ScorePair,DeNovoRecord> motifPairsMap, 
			Set<Gene> expGenesCutoff0, Map<Gene, Set<RecS>> discTableA){
		
		for (Map.Entry<ScorePair,DeNovoRecord> p : motifPairsMap.entrySet()){
			
			String motif = p.getKey().getKey();
			Double score = p.getKey().getValue();
			
			
			
			int bls = p.getValue().getBlsMin();
			int F = p.getValue().getNumFamilies();
			double conf = p.getValue().getConfidence();
			
			Set<Gene> motifTargets = getMotifTargetsFromTable(discTableA,motif);
			
			//Set<Gene> intersBoldA = calculateIntersection(motifTargets,expGenesNoCutoff);
			Set<Gene> intersBold = calculateIntersection(motifTargets,expGenesCutoff0);
			
			StringBuilder entry = new StringBuilder();
			entry.append(motif + "\t");
			entry.append(bls+"\t");
			entry.append(F+"\t");
			entry.append(conf+"\t");
			entry.append(score + "\t");
			entry.append(motifTargets.size()+"\t");
			//entry.append(intersBoldA.size()+"\t");
			entry.append(intersBold.size()+"\t");
			entry.append(intersBold);
			System.out.println(entry.toString());

			

		}
		
	}
	
	private Set<Gene> getMotifTargetsFromTable(Map<Gene, Set<RecS>> table,
			String motif) {
		
		Set<Gene> genes = new HashSet<Gene>();
		
		
		for (Map.Entry<Gene, Set<RecS>> e : table.entrySet()){
		
			inner:for (RecS record : e.getValue()){
				
			
				String m = record.getBindingSite();
				
				if (!m.equals(motif)){ //skip all records not matching with the motif
					continue;
				}
				
				//matches with motif
				//System.out.println("match");
				genes.add(e.getKey());
				break inner;
				
			}
			
		}
		return genes;
		
	}

	private SortedMap<ScorePair, DeNovoRecord> generateMotifMapFromDeNovoRecords(
			Set<DeNovoRecord> kn1Variants) {
		
		SortedMap<ScorePair, DeNovoRecord> map = new TreeMap<ScorePair, DeNovoRecord>();
		
		for (DeNovoRecord record : kn1Variants){
						
			ScorePair p = new ScorePair(record.getMotif(),record.getPwmScore());
			map.put(p,record);
			
		}
		
		return map;
		
	}

	/**
	 * Different Collections:
	 * Bolduc and Bolduc_O , Bolduc_O C Bolduc!
	 * PWMScan
	 * Discovery
	 * @param pwmThreshold
	 * @throws IOException 
	 */
	public void generateVennDiagram(String exp, int t) throws IOException{
		
		
		System.out.println("Venndiagram for: "+t);
		System.out.println("");
	
		Map<Gene, Set<RecS>> bolducTable = KN1Toolbox.readTableSFromFile(KN1Toolbox.generateOutputFilenameBolduc(t));
		Map<Gene, Set<RecF>> bolducOTable = KN1Toolbox.readTableFFromFile(KN1Toolbox.generateOutputFilenameBolducWithOrthologs(t));
		Map<Gene, Set<RecS>> pwmScTable = KN1Toolbox.readTableSFromFile(KN1Toolbox.generateOutputFilenamePWMScan(t));
		
		Map<Gene, Set<RecS>> discTable;
		discTable = KN1Toolbox.readTableSFromFile(KN1Toolbox.generateOutputFilenameDeNovoWorstBS(exp,t));

		//all intersections:
		
		//C1_4
		Set<Gene> Bolduc = new HashSet<Gene>(bolducTable.keySet());
		Set<Gene> BolducO = new HashSet<Gene>(bolducOTable.keySet());
		Set<Gene> PWMSC = new HashSet<Gene>(pwmScTable.keySet());
		Set<Gene> Disc = new HashSet<Gene>(discTable.keySet());
		
		System.out.println("Bolduc: "+Bolduc.size());
		System.out.println("BolducO: "+BolducO.size());
		System.out.println("PWMSC: "+PWMSC.size());
		System.out.println("Disc: "+Disc.size());
		
		/*System.out.println("Bolduc: "+Bolduc);
		System.out.println("BolducO: "+BolducO);
		System.out.println("Bolduc: "+BolducO);
		System.out.println("Bolduc: "+BolducO);*/

		

		
		//C2_4
		Set<Gene> Bolduc_I_BolducO = calculateIntersection(Bolduc,BolducO); //12
		Set<Gene> Bolduc_I_PWMSC = calculateIntersection(Bolduc,PWMSC);//13
		Set<Gene> Bolduc_I_Disc = calculateIntersection(Bolduc,Disc);//14
		Set<Gene> BolducO_I_PWMSC = calculateIntersection(BolducO,PWMSC);//23
		Set<Gene> BolducO_I_Disc = calculateIntersection(BolducO,Disc);//24
		Set<Gene> PWMSC_I_Disc = calculateIntersection(PWMSC,Disc);//34
		
		System.out.println("Bolduc_I_BolducO: "+ Bolduc_I_BolducO.size());
		System.out.println("Bolduc_I_PWMSC: "+ Bolduc_I_PWMSC.size());
		System.out.println("Bolduc_I_Disc: "+ Bolduc_I_Disc.size());
		System.out.println("BolducO_I_PWMSC: "+ BolducO_I_PWMSC.size());
		System.out.println("BolducO_I_Disc: "+ BolducO_I_Disc.size());
		System.out.println("PWMSC_I_Disc: "+ PWMSC_I_Disc.size());
		
		
		//C3_4
		Set<Gene> Bolduc_I_BolducO_I_PWMSC = calculateIntersection(Bolduc_I_BolducO, PWMSC); //123
		Set<Gene> Bolduc_I_BolducO_I_Disc = calculateIntersection(Bolduc_I_BolducO, Disc); //124
		Set<Gene> Bolduc_I_PWMSC_I_Disc = calculateIntersection(Bolduc_I_PWMSC, Disc); //134
		Set<Gene> BolducO_I_PWMSC_I_Disc = calculateIntersection(BolducO_I_PWMSC, Disc); //234
		
		System.out.println("Bolduc_I_BolducO_I_PWMSC: "+ Bolduc_I_BolducO_I_PWMSC.size());
		System.out.println("Bolduc_I_BolducO_I_Disc: "+ Bolduc_I_BolducO_I_Disc.size());
		System.out.println("Bolduc_I_PWMSC_I_Disc: "+ Bolduc_I_PWMSC_I_Disc.size());
		System.out.println("BolducO_I_PWMSC_I_Disc: "+ BolducO_I_PWMSC_I_Disc.size());
		
		
		//C4_4
		Set<Gene> FullIntersection = calculateIntersection(Bolduc_I_BolducO, PWMSC_I_Disc);
		System.out.println("FullIntersection: "+ FullIntersection.size());
		
		System.out.println("");
		
		System.out.println(FullIntersection);
		
		
		System.out.println("DiffBolducOenFullInters: "+calculateDifference(FullIntersection,BolducO));
		
		
		System.out.println("Extra discovery targets are FP or functional?");
		System.out.println("Intersection: "+FullIntersection);
		System.out.println("");
		System.out.println("SDISC only: "+calculateDifference(FullIntersection, Disc));
		System.out.println("");
		System.out.println("Bolduc zonder BolducO: "+calculateDifference(BolducO,Bolduc));

		
		
	}

	private Set<Gene> calculateIntersection(Set<Gene> s1, Set<Gene> s2) {
		Set<Gene> intersection = new HashSet<Gene>();
		for (Gene g : s1){
			if (s2.contains(g)){
				intersection.add(g);
			}
		}
		return intersection;
		
	}
	
	private Set<Gene> calculateDifference(Set<Gene> s1, Set<Gene> s2) {
		Set<Gene> diff = new HashSet<Gene>();
		for (Gene g : s1){
			if (!s2.contains(g)){
				diff.add(g);
			}
		}
		for (Gene g : s2){
			if (!s1.contains(g)){
				diff.add(g);
			}
		}
		return diff;
		
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
	
		//KN1Toolbox.setKN1Directory("/home/ddewitte/Bureaublad/CloudSpellerExperimentsFinal/ABAnalysis/");

		
		String [] experiments = {"C50F2","C50F3","C50F4","C50F5","C50F10","C50F15","C50F20","C50F25","C50F50",
				"C80F2","C80F3","C80F4","C80F5","C80F10","C80F15","C80F20","C80F25","C80F50",
				"C90F2","C90F3","C90F4","C90F5","C90F10","C90F15","C90F20","C90F25","C90F50"
		};
		
		for (int i=0; i<experiments.length; i++){
		
			KN1Toolbox.setExperiment(experiments[i]);
			
			KN1OverlapAnalyzer overlapAnalyzer = new KN1OverlapAnalyzer();
	
			System.out.println("exp= "+KN1Toolbox.getExperiment());
			overlapAnalyzer.generateVennDiagram(KN1Toolbox.getExperiment(),KN1Toolbox.getPWMScoreThreshold());
			
			overlapAnalyzer.singleMotifOverlap(KN1Toolbox.getExperiment());
		}	

	}





}
