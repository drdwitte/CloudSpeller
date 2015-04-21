package be.iminds.cloudspeller.IPA1Analysis;

import be.iminds.cloudspeller.input.Gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;


import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifpermutationgroups.IUPACContent;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;
import be.iminds.cloudspeller.motifpermutationgroups.MotifPermutationGroup;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;

import be.iminds.cloudspeller.KN1Analysis.PMRecord;
import be.iminds.cloudspeller.KN1Analysis.RecS;
import be.iminds.cloudspeller.TargetMatching.PermutationPatternFileGenerator;

public class DeNovoMotifMatching {

	private static String [] motifs;
	static {
		int [] thresholds = {15/*,50,60,70,90*/,95};
		BLS.initializeBLSConstants(thresholds);
		FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
		 motifs = new String[3];
		 motifs[0]="TGGGCC";
		 motifs[1]="TGGGCT";
		 motifs[2]="TGGGCY";

	}
	
	public static void generateTMFiles() throws IOException{
		
		
		
		for (String motif : motifs){
		
			for (int bls : BLS.getBLSThresholds()){
				PermutationPatternFileGenerator.generateTMFile(motif,bls);
			}
		}
		
	}
	
	public static void main (String [] args) throws IOException{
		
		//generateTMFiles();
		//reworkPMOutputs("AF");
		//reworkPMOutputs("AB");
		
		reworkOCRPMs("AB");
		reworkOCRPMs("AF");
		
		//generateConfidenceChartsAndDistributions("AB");
		//generateConfidenceChartsAndDistributions("AF");
		
		generateGeneMotifTables("AF");
		generateGeneMotifTables("AB");
		
		
		
		
	}

	
	
	private static void reworkOCRPMs(String type) throws IOException{

		System.out.println("rework "+type);
		
		
		for (String motif : motifs){
			for (int bls : BLS.getBLSThresholds()){
				
				String filename = IPA1Toolbox.getIPA1Directory()+"IPA1OCRFiltering/OCR_filtered_PMs"
							+type+motif+bls+".bed";					
						
				System.out.println(type+motif+bls);
				BufferedReader reader = new BufferedReader(new FileReader(filename));
				String outputFilename=IPA1Toolbox.getPMFilename(type, motif, bls);	
				BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilename));
				String line;
				
				while ( (line = reader.readLine())!=null){		
				
					String [] spl = line.split("\t");
					
					if (spl.length!=9)System.out.println(line);
					String newLine="";
					for (int i=3; i<spl.length-1; i++){
						newLine+=spl[i]+"\t";
					}
					newLine+=spl[spl.length-1];
					
					writer.append(newLine); writer.newLine();
					
				}
				
				writer.close();
				
			}
		}
		
		
		
		
		
		
	}

	@SuppressWarnings("unused")
	private static void reworkPMOutputs(String type) throws IOException {
		System.out.println("rework "+type);
		String filename = IPA1Toolbox.getIPA1Directory()+type+"PM.txt";
		
		
		int c=0;
		for (String motif : motifs){
			
			for (int bls  : BLS.getBLSThresholds()){
				System.out.println(type+motif+bls);
				BufferedReader reader = new BufferedReader(new FileReader(filename));
				String outputFilename=IPA1Toolbox.getPMFilename(type, motif, bls);	
				BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilename));
				String line;
				
				while ( (line = reader.readLine())!=null){
					
					
					//record.set(line);
					c++;
					if (c%1000000==0){
						if (c%10000000==0){
							System.out.println("*");
						} else {
							System.out.print("*");
						}
					}
					if (line.startsWith(motif+"_"+bls)){
					//if (record.getBls()==bls && record.getMotif().equals(motif)){
						writer.write(line); writer.newLine();
					}
				}
				reader.close();
				writer.close();
			}
		}
		
			
	}
	
	@SuppressWarnings("unused")
	private static void reworkPMOutputsFull(String type) throws IOException {
		System.out.println("rework "+type);
		String filename = IPA1Toolbox.getIPA1Directory()+type+"PMFull.txt";
		
		
		int c=0;
		for (String motif : motifs){
			
			for (int bls  : BLS.getBLSThresholds()){
				System.out.println(type+motif+bls);
				BufferedReader reader = new BufferedReader(new FileReader(filename));
				String outputFilename=IPA1Toolbox.getPMFilename(type, motif, bls);	
				BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilename));
				String line;
				
				while ( (line = reader.readLine())!=null){
					
					
					//record.set(line);
					c++;
					if (c%1000000==0){
						if (c%10000000==0){
							System.out.println("*");
						} else {
							System.out.print("*");
						}
					}
					if (line.startsWith(motif+"_"+bls)){
					//if (record.getBls()==bls && record.getMotif().equals(motif)){
						writer.write(line); writer.newLine();
					}
				}
				reader.close();
				writer.close();
			}
		}
		
			
	}
	
	
	@SuppressWarnings("unused")
	private static void generateConfidenceChartsAndDistributions(String type) throws IOException {
		System.out.println("confCharts "+type);
		String filename = IPA1Toolbox.getIPA1Directory()+type+"PM.txt";
		
		Map<String,Integer> motifBLSFreqMap = new HashMap<String,Integer>();
		
		int c=0;

		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String line;
			
		Map<String,Set<String>> motifBLSFamMap = new HashMap<String,Set<String>>();
		
		while ( (line = reader.readLine())!=null){
			
			c++;
			if (c%1000000==0){
				if (c%10000000==0){
					System.out.println("*");
				} else {
					System.out.print("*");
				}
			}
			
			Scanner scanner = GeneralToolbox.generateScanner(line);
			
			String motifBLS = scanner.next(); //motifBLS
			
			String fam = scanner.next(); //family
			
			Set<String> fams = motifBLSFamMap.get(motifBLS);
			
			if (fams==null){
				fams = new HashSet<String>();
				motifBLSFamMap.put(motifBLS,fams);
			}
			
			fams.add(fam);
		}
			
		//motifBLSFamMap.clear();
		reader.close();
		
		System.out.println("");
		
		for (Map.Entry<String,Set<String>> e : motifBLSFamMap.entrySet()){
			motifBLSFreqMap.put(e.getKey(),e.getValue().size());
		}
		motifBLSFamMap.clear();
		
		
		emitConfCharts(motifBLSFreqMap);
		generateDistribution(motifBLSFreqMap);
			
			
			
			
			
			
			
		


			
	}

	private static void generateDistribution(
			Map<String, Integer> motifBLSFreqMap) {
		
		//int interv = 50;
		
		Map<String,ArrayList<Integer>> contentOccs = new HashMap<String,ArrayList<Integer>>();
		IUPACContent content = new IUPACContent("test");
		
		for (Map.Entry<String,Integer> e : motifBLSFreqMap.entrySet()){
			
			String motifBLS = e.getKey();
			String [] spl = motifBLS.split("_");
			content.setContentFromStr(spl[0]);
			String key = content.toString()+"_"+spl[1];
			
			ArrayList<Integer> oldValue = contentOccs.get(key);
			
			if (oldValue == null){
				oldValue = new ArrayList<Integer>();
				contentOccs.put(key,oldValue);
				
			}
			
			oldValue.add(e.getValue());
			
		}
		
		

		for (Map.Entry<String,ArrayList<Integer>> e : contentOccs.entrySet()){

			System.out.println(e.getKey());
			
			if (e.getKey().contains("95")){
				
			} else {
				continue;
			}
			
			System.out.println(e.getValue());
			
			/*SortedMap<Integer,Integer> histogram = new TreeMap<Integer,Integer>();


			for (Integer i : e.getValue()){
				
				int ID = i/interv;
				
				Integer freq = histogram.get(ID);
				
				if (freq == null){
					histogram.put(ID,1);
				} else {
					histogram.put(ID,freq+1);
				}
				
			}
			
			for (Map.Entry<Integer,Integer> h : histogram.entrySet()){
				int ID = h.getKey();
				int start = ID*interv;
				System.out.println((start+interv/2)+"\t"+h.getValue());
				
			}*/
			
			

		}
		
		for (Map.Entry<String,Integer> e : motifBLSFreqMap.entrySet()){
			System.out.println(e.getKey()+"\t"+e.getValue());
		}

		System.out.println("");
		System.out.println("");
		System.out.println("");


		
		
	}

	private static void emitConfCharts(Map<String,Integer> motifBLSFreqMap) {
		//conf charts for ref motifs
		for (String motif : motifs){
			
			
			MotifContent content = new IUPACContent(motif);
			Set<Motif> permutations = content.createPermutationGroup(1000);

			
			System.out.println(motif);
			
			
			
			for (int bls : BLS.getBLSThresholds()){
			
				ArrayList<Integer> occs = new ArrayList<Integer>();
				
				Integer F = motifBLSFreqMap.get(motif+"_"+bls);
				
				for (Motif p : permutations){
					
					Integer o = motifBLSFreqMap.get(p.toString()+"_"+bls);
					occs.add((o==null)?0:o);
				}
				
				double Fb = MotifPermutationGroup.findMedian(occs);

				double C = 100.0 * (F-Fb)/F;
				System.out.println(bls+"\t"+F+"\t"+Fb+"\t"+C);
				
				
			}
		}
			
			
	}

	private static void generateGeneMotifTables(String type) throws IOException {
		
		for (int bls : BLS.getBLSThresholds()){
		
			Map<Gene,Set<RecS>> table = new HashMap<Gene,Set<RecS>>(); 
			for (String motif : motifs){	
				String filename = IPA1Toolbox.getPMFilename(type,motif,bls);
				BufferedReader reader = new BufferedReader(new FileReader(filename));
				
				String line;
				//AACGTCDACGKG_15	iORTHO033607	SB01G017030	-	-220	-209
				Set<PMRecord> pmRecords = new HashSet<PMRecord>();
				while ( (line=reader.readLine()) !=null)
				{
					//DEBUG extraatje:
					/*String [] spl = line.split("\t");
					int start=Integer.parseInt(spl[spl.length-2]);
					if (start>=-500)*/
					System.out.println(line);
					pmRecords.add(new PMRecord(line));
				}
				reader.close();
				
				for (PMRecord record : pmRecords){
					if (record.getGene().startsWith(IPA1Toolbox.species)){
						
						String g = record.getGene();
						Gene gene = new Gene(g,IPA1Toolbox.species);
						
						Set<RecS> set = table.get(gene);
						
						if (set==null){
							set = new HashSet<RecS>();
							table.put(gene,set);
						}
						
						set.add(new RecS(record.getMotif(),-1));
						
						
					}
				}
				
				
				
				
			}
			IPA1Toolbox.printTableS(table,IPA1Toolbox.getDeNovoFilename(type,bls));
		}
		
		
	}
	


	
	
}
