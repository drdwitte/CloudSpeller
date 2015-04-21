package be.iminds.cloudspeller.TargetMatching;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.toolbox.LineIterator;

import be.iminds.cloudspeller.KN1Analysis.PMRecord;

import be.iminds.cloudspeller.motifpermutationgroups.MotifPermutationGroup;



public class PMOutputToConfCharts {

	public static void main (String [] args) throws IOException{
		
		//FIXME moet sampelen om conf charts te kunnen bepalen want kan ook 0 x voorkomen!!
		System.out.println("WArning fixme!!");
		
		String outputPM = args[0];
		String confOutput = args[1];
		String queryMotif = args[2];
		
		//ACATGNATG_40	iORTHO015554	BD5G10200	-	-1559	-1551
		LineIterator it = new LineIterator(outputPM);
		
		Map<String,Set<String>> motifFamMap = new HashMap<String,Set<String>>();
		
		
		while (it.hasNext()){
			String line = it.next();
			PMRecord record = new PMRecord(line);
			
			String motif = record.getMotif();
			Set<String> famMatches = motifFamMap.get(motif);
			
			if (famMatches == null){
				famMatches = new HashSet<String>();
				motifFamMap.put(motif,famMatches);
			}
			famMatches.add(record.getFamily());
			
			
		}

		ArrayList<Integer> familyOccs = new ArrayList<Integer>();

		
		for (Map.Entry<String,Set<String>> e : motifFamMap.entrySet()){
			//System.out.println(e.getKey()+"\t"+e.getValue().size() + "\t #degPos= "+getNumDegPos(e.getKey()));
			familyOccs.add(e.getValue().size()); 
		}
		
		

		BufferedWriter writer = new BufferedWriter(new FileWriter(confOutput));
		double median = MotifPermutationGroup.findMedian(familyOccs);
		writer.append("Median is: "+median);
		writer.newLine();
		
		if (motifFamMap.get(queryMotif)==null){
			
			writer.append("Query motif: "+queryMotif+"\t NOT FOUND");
			
		} else {
		
			double F = motifFamMap.get(queryMotif).size();
			double C = (F - median)/F; 
			writer.append("Query motif: "+queryMotif+"\t"+F+"\t"+C);
		}
		
		
		writer.newLine();
		writer.newLine();
		for (Map.Entry<String,Set<String>> e : motifFamMap.entrySet()){
			double F = e.getValue().size();
			double C = (F - median)/F; 
			
			writer.append(e.getKey()+"\t"+F+"\t"+C);
			writer.newLine();
			/*for (String s : e.getValue()){
				System.out.print(s+" ");
			}
			System.out.println("");
			System.out.println("");*/
		}
		
		writer.close();
		
		
	}
	
	
	@SuppressWarnings("unused")
	private static int getNumDegPos(String s){
		int n = 0;
		for (int i=0; i<s.length(); i++){
			char c = s.charAt(i);
			
			if (c=='A' || c=='C' || c=='G' || c=='T'){
				continue;
			}
			n++;
		}
		return n;
	}
}
