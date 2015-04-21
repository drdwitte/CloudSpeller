package be.iminds.cloudspeller.KN1Analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;




import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import be.iminds.cloudspeller.toolbox.LineIterator;


public class RedundancyGroupClustering {

	public static void main (String [] args) throws IOException{
		String expDir = "/home/ddewitte/Bureaublad/CloudSpellerExperimentsFinal/KN1_6NOV/";
		String file = "RedundancyClustering/RedundancyC90F0.txt";
		
		
		
		Map<String,Set<Variant>> geneMatchingVariantsMap  = importVariants(expDir+file);
		
		
		
		int minF=0;
		
		
		Set<Variant> essentialVariants = reduceVariants(geneMatchingVariantsMap,minF);

		System.out.println("Essential variants: "+essentialVariants);
		
		
		StringBuilder [] charsAtPos = {new StringBuilder(), new StringBuilder(), new StringBuilder(), 
				new StringBuilder(), new StringBuilder(), new StringBuilder(), new StringBuilder(), 
				new StringBuilder(), new StringBuilder(), new StringBuilder(), new StringBuilder(), 
				new StringBuilder() };
		
		
		for (Variant v : essentialVariants){
			for (int i=0; i<v.length(); i++){
				charsAtPos[i].append(v.charAt(i));
			}
		}
		
		for (int i=0; i<charsAtPos.length; i++){
			System.out.println(charsAtPos[i]);
		}
		
		createVisualization(essentialVariants, geneMatchingVariantsMap);
		

		
		
			
			
			
			
		
		
		
	}
	
	private static void createVisualization(Set<Variant> essentialVariants,
			Map<String, Set<Variant>> geneMatchingVariantsMap) {
		
		Set<Variant> allVariants = new HashSet<Variant>();
		
		Map<String,Integer> gIDMap = new HashMap<String, Integer>();
		
		
		SortedMap<Integer,Set<String>> numberOfMatchesPerGeneMap = new TreeMap<Integer,Set<String>>(java.util.Collections.reverseOrder());
		
		//sort genes according to number of variant matches to improve cluster visualization
		for (Map.Entry<String,Set<Variant>> e : geneMatchingVariantsMap.entrySet()){
			allVariants.addAll(e.getValue());
			Set<String> genes = numberOfMatchesPerGeneMap.get(e.getValue().size());
			if (genes == null){
				genes = new HashSet<String>();
				numberOfMatchesPerGeneMap.put(e.getValue().size(),genes);
			}
			genes.add(e.getKey());
		}
		
		int gId=0;
		
		System.out.println("GeneIDS:");
		for (Map.Entry<Integer,Set<String>> e : numberOfMatchesPerGeneMap.entrySet()){
			for (String s : e.getValue()){
				gIDMap.put(s,gId++);
				System.out.println(s+"\t"+(gId-1));
			}
		}
		
		
		
		System.out.println("");
		System.out.println("");
		System.out.println("All Variants");
		System.out.println("");
		Map<Variant,Integer> vIDMap = new HashMap<Variant, Integer>();
		int vID=0;
		for (Map.Entry<String,Set<Variant>> e : geneMatchingVariantsMap.entrySet()){
			
			for (Variant v : e.getValue()){
				if (allVariants.contains(v)){
					vIDMap.put(v,vID++);
					allVariants.remove(v);
					System.out.println(v+"\t"+(vID-1));
				}
			}
		}
				
		boolean [][] table = new boolean [vIDMap.size()][gIDMap.size()];
		
		
		for (Map.Entry<String,Set<Variant>> e : geneMatchingVariantsMap.entrySet()){
			for (Variant v : e.getValue()){
			
				int vI=vIDMap.get(v);
				int gI=gIDMap.get(e.getKey());
				table[vI][gI]=true;
			}
		}

		System.out.println("");
		System.out.println("Essential variants: ");
		System.out.println("");
		ArrayList<Integer> essIDs = new ArrayList<Integer>();
		for (Variant v : essentialVariants){
			System.out.println(v + "\t" + vIDMap.get(v));
			essIDs.add(vIDMap.get(v));
		}
		
		
		
		
		
		
		System.out.println("");
		System.out.println("Clustering Matrix");
		
		
		StringBuilder x = new StringBuilder();
		StringBuilder y = new StringBuilder();
		StringBuilder xess = new StringBuilder();
		StringBuilder yess = new StringBuilder();
		
		
		
		
		for (int i=0; i<table.length; i++){
			
			for (int j=0; j<table[i].length; j++){
								
				if (table[i][j]){
					if (essIDs.contains(i)){
						xess.append(i+"\n");
						yess.append(j+"\n");
					} else {
						x.append(i+"\n");
						y.append(j+"\n");
					}
				}
			}
		}
		
		
		System.out.println("Xcoordinaten");
		System.out.println(x);
		
		System.out.println("Ycoordinaten");
		System.out.println(y);
		
		System.out.println("XcoordinatenEss");
		System.out.println(xess);
		
		System.out.println("YcoordinatenEss");
		System.out.println(yess);
		
		
		
		

		
		


		
		
		
		
	}

	private static Set<Variant> reduceVariants(Map<String,Set<Variant>> geneMatchingVariantsMap,int minF) {

		Set<Variant> essentialVariants = new HashSet<Variant>();
		
		for (Map.Entry<String,Set<Variant>> e : geneMatchingVariantsMap.entrySet()){
			
			Iterator<Variant> varIt = e.getValue().iterator();
			Variant bestVar = varIt.next();
			
			
			while (varIt.hasNext()){
				Variant v = varIt.next();
				if (v.getF()>bestVar.getF()){
					bestVar=v;
				}
				
			}
			
			if (bestVar.getF()>=minF){
				essentialVariants.add(bestVar);
			}
			
			
			
		}
		
		int numBolducOTargets=0;
		
		for (Map.Entry<String,Set<Variant>> e : geneMatchingVariantsMap.entrySet()){
			int numVars=0;
			for (Variant var : e.getValue()){
			
				if (essentialVariants.contains(var)){
					System.out.println(e.getKey()+"\t"+var);
					numVars++;
				}
				
			}
			
			if (numVars>0){
				numBolducOTargets++;
			}
			
			
			System.out.println("");

			
		}
		
		System.out.println("Number of bolduc targets with F>="+minF+" : "+numBolducOTargets);
		System.out.println("Number of essential variants: "+essentialVariants.size());
		
		return essentialVariants;
		
	}

	public static Map<String,Set<Variant>> importVariants(String variantsFile) throws IOException{
		LineIterator iterator = new LineIterator(variantsFile);

		Map<String,Set<Variant>> geneMatchingVariantsMap = new HashMap<String,Set<Variant>>();
		
		
		while (iterator.hasNext()){
			
			//variant 	 blsMin 	 #Fam 	 Conf 	 score 	 #maize targets 	  [ovBold]
			
			String line = iterator.next();
			//System.out.println(line);
			Scanner scanner = GeneralToolbox.generateScanner(line);
			String variantStr = scanner.next();
			scanner.next();
			int F = Integer.parseInt(scanner.next());
			Variant variant = new Variant(variantStr,F);
			
			String targetsStr = line.substring(line.indexOf('[')+1,line.indexOf(']'));
			
			if (targetsStr.length()==0){
				continue;
			}
			
			String [] splits = targetsStr.split(",");
			
			
			
			for (String target : splits){
				
				
				addTargetVariant(target.trim(),variant,geneMatchingVariantsMap);
			}
		}
		
		
		return geneMatchingVariantsMap;
	}
	

	private static void addTargetVariant(String target, Variant variant, Map<String,Set<Variant>> geneMatchingVariantsMap) {
		Set<Variant> variants = geneMatchingVariantsMap.get(target);
		if (variants==null){
			variants = new HashSet<Variant>();
			geneMatchingVariantsMap.put(target,variants);
		}
		variants.add(variant);
	}
}


class Variant {

	private String motif;
	private int F;
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Variant){
			return this.toString().equals(obj.toString());
		}
		return false;
		
	}
	
	public int length() {
		return motif.length();
	}

	public char charAt(int i) {
		return motif.charAt(i);
	}

	public int getF(){
		return F;
	}


	@Override
	public int hashCode() {
		return this.toString().hashCode();
	}

	@Override
	public String toString() {
		return motif+F;
	}


	
	public Variant(String motif, int F){
		this.motif=motif;
		this.F=F;
		
	}
	
	
	
	
}