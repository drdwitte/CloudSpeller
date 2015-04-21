package be.iminds.cloudspeller.IPA1Analysis;

import be.iminds.cloudspeller.input.Gene;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.FreqVec;

import be.iminds.cloudspeller.phylogenetics.BLS;

import be.iminds.cloudspeller.KN1Analysis.RecF;
import be.iminds.cloudspeller.KN1Analysis.RecS;


public class OverlapAnalyzer {

	private static String [] motifs;
	static {
		int [] thresholds = {15,95};
		BLS.initializeBLSConstants(thresholds);
		FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
		 motifs = new String[3];
		 motifs[0]="TGGGCC";
		 motifs[1]="TGGGCT";
		 motifs[2]="TGGGCY";

	}
	
	public static void main (String [] args) throws IOException{
		
		String [] types = {"AF","AB"};
		for (String type : types){
			for (int bls : BLS.getBLSThresholds()){
				generateVennDiagram(type,bls);
				
				for (String motif : motifs){
					
					singleMotifOverlap(type,motif, bls);
				}
				
			}
		}
	}
	
	
	public static void generateVennDiagram(String type, int bls) throws IOException{
		
		System.out.println("Venndiagram for Type "+type+"BLS Ti= "+bls);
		System.out.println("");
	
		Map<Gene, Set<RecS>> expTable = IPA1Toolbox.readTableSFromFile(IPA1Toolbox.getFilenameExperimentalMatches());
		Map<Gene, Set<RecF>> expOTable = IPA1Toolbox.readTableFFromFile(IPA1Toolbox.getFilenameExperimentalMatchesWithOrthologs());
		Map<Gene, Set<RecS>> pwmScTable = IPA1Toolbox.readTableSFromFile(IPA1Toolbox.getFilenamePWMMatches());
		
		Map<Gene, Set<RecS>> discTable;
		discTable = IPA1Toolbox.readTableSFromFile(IPA1Toolbox.getDeNovoFilename(type,bls));

		//all intersections:
		
		//C1_4
		Set<Gene> Exp = new HashSet<Gene>(expTable.keySet());
		Set<Gene> ExpO = new HashSet<Gene>(expOTable.keySet());
		Set<Gene> PWMSC = new HashSet<Gene>(pwmScTable.keySet());
		Set<Gene> Disc = new HashSet<Gene>(discTable.keySet());
		
		System.out.println("Exp: "+Exp.size());
		System.out.println("ExpO: "+ExpO.size());
		System.out.println("PWMSC: "+PWMSC.size());
		System.out.println("Disc: "+Disc.size());
		
		//C2_4
		Set<Gene> Exp_I_ExpO = IPA1Toolbox.calculateIntersection(Exp,ExpO); //12
		Set<Gene> Exp_I_PWMSC = IPA1Toolbox.calculateIntersection(Exp,PWMSC);//13
		Set<Gene> Exp_I_Disc = IPA1Toolbox.calculateIntersection(Exp,Disc);//14
		Set<Gene> ExpO_I_PWMSC = IPA1Toolbox.calculateIntersection(ExpO,PWMSC);//23
		Set<Gene> ExpO_I_Disc = IPA1Toolbox.calculateIntersection(ExpO,Disc);//24
		Set<Gene> PWMSC_I_Disc = IPA1Toolbox.calculateIntersection(PWMSC,Disc);//34
		
		System.out.println("Exp_I_ExpO: "+ Exp_I_ExpO.size());
		System.out.println("Exp_I_PWMSC: "+ Exp_I_PWMSC.size());
		System.out.println("Exp_I_Disc: "+ Exp_I_Disc.size());
		System.out.println("ExpO_I_PWMSC: "+ ExpO_I_PWMSC.size());
		System.out.println("ExpO_I_Disc: "+ ExpO_I_Disc.size());
		System.out.println("PWMSC_I_Disc: "+ PWMSC_I_Disc.size());
				
		//C3_4
		Set<Gene> Exp_I_ExpO_I_PWMSC = IPA1Toolbox.calculateIntersection(Exp_I_ExpO, PWMSC); //123
		Set<Gene> Exp_I_ExpO_I_Disc = IPA1Toolbox.calculateIntersection(Exp_I_ExpO, Disc); //124
		Set<Gene> Exp_I_PWMSC_I_Disc = IPA1Toolbox.calculateIntersection(Exp_I_PWMSC, Disc); //134
		Set<Gene> ExpO_I_PWMSC_I_Disc = IPA1Toolbox.calculateIntersection(ExpO_I_PWMSC, Disc); //234
		
		System.out.println("Exp_I_ExpO_I_PWMSC: "+ Exp_I_ExpO_I_PWMSC.size());
		System.out.println("Exp_I_ExpO_I_Disc: "+ Exp_I_ExpO_I_Disc.size());
		System.out.println("Exp_I_PWMSC_I_Disc: "+ Exp_I_PWMSC_I_Disc.size());
		System.out.println("ExpO_I_PWMSC_I_Disc: "+ ExpO_I_PWMSC_I_Disc.size());
		
		
		//C4_4
		Set<Gene> FullIntersection = IPA1Toolbox.calculateIntersection(Exp_I_ExpO, PWMSC_I_Disc);
		System.out.println("FullIntersection: "+ FullIntersection.size());
		

	}

private static void singleMotifOverlap(String type,String motif, int bls) throws IOException {

	Map<Gene, Set<RecF>> expOTable = IPA1Toolbox.readTableFFromFile(IPA1Toolbox.getFilenameExperimentalMatchesWithOrthologs());
	Map<Gene, Set<RecS>> discTable;
	discTable = IPA1Toolbox.readTableSFromFile(IPA1Toolbox.getDeNovoFilename(type,bls));
	//System.out.println(discTable.size());
	removeOtherMotifs(motif,discTable);
	//System.out.println(discTable.size());
	Set<Gene> intersection = IPA1Toolbox.calculateIntersection(discTable.keySet(),expOTable.keySet());
		
	StringBuilder entry = new StringBuilder();
	entry.append(motif + "\t");
	entry.append(bls+"\t");
	entry.append(discTable.size()+"\t");
	entry.append(intersection.size()+"\t");
	//entry.append(intersection);
	System.out.println(entry.toString());
	
	
	
	
}

private static  void removeOtherMotifs(String motif, Map<Gene, Set<RecS>> discTable) {

for (Map.Entry<Gene,Set<RecS>> e : discTable.entrySet()){
		
		Set<RecS> filtered = new HashSet<RecS>();
		for (RecS rec : e.getValue()){
			if (rec.getBindingSite().equals(motif)){
				filtered.add(rec);
			}
		}
		
		e.setValue(filtered);
	}
	
	Set<Gene> keysToRemove = new HashSet<Gene>();
	for (Map.Entry<Gene,Set<RecS>> e : discTable.entrySet()){
		if (e.getValue().size()==0){
			keysToRemove.add(e.getKey());
		}
	}
	
	for (Gene g : keysToRemove){
		discTable.remove(g);
	}
	
}


}
