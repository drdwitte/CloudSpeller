package be.iminds.cloudspeller.KN1Analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import be.iminds.cloudspeller.toolbox.LineIterator;

public class PostProcessMaizeGeneAnalysis {

	private static final String Bolduc = "B";
	private static final String OrthoMatch = "O"; 
	private static final String PWMSCMaize = "SM";
	private static final String PWMSCOrtho = "SO"; 
	

	
	
	public static void main (String [] args) throws IOException{
		
		
		
		//class - Gene
		Map<String,Set<String>> classGeneMap = initializeClassGeneMap();
		
		String filename = "/home/ddewitte/Bureaublad/CloudSpellerExperiments/NieuweKN1Analyse/outputDetailedKN1Analysis.txt";
		
		LineIterator it = new LineIterator(filename);
		
		while (it.hasNext()){
			//read entry
			String line = it.next();
			if (line.length()==0 || !line.equals("*")){
				continue;
			}
				
			line = it.next();
			Scanner s = createScanner(line);
			
			String gene = s.next();
			boolean b = s.next().charAt(0) == 'y';
			s.close();
			
			line=it.next(); //BS in maize
			
					
			
			while (!line.startsWith("Max")){
				line = it.next();
				
			}
			
			s = createScanner(line);
			s.next();
			boolean pm = s.nextDouble() >= 0.0; //check maize pwm score
			s.close();
			
			line = it.next(); //get genefamily line
			
			
			line = it.next();
			while (!line.startsWith("Highest")){ //skip ortholog binding sites
				line = it.next();
			}
			
			s = createScanner(line);
			s.next();
			double orthoScore = s.nextDouble();
			boolean o = orthoScore > (MaisGeneAnalyzer.lowerBoundScore+1);
			boolean po = orthoScore >=0.0;
			s.close();
			
			classGeneMap.get(generateKeyForClass(b, o, pm, po)).add(gene);
			
			
			//System.out.println(generateKeyForClass(b, o, pm, po)+gene);
			
			//it.next(); //Least Degenerate...
			//it.next(); //motifs
			
			/*
			ZM03G25960	n
			BS in maize Gene:
			TGAGTGACCGAC	-26.409208775416467
			Max maize score: -26.409208775416467
			GeneFamilies: 1
			iORTHO004412
			SB03G040860	TGAGAGATGGAG	-28.917612580390063
			BD2G55980	TGATCGATGGAG	-10.132288722849788
			OS01G64590	TGATTGATGGAG	-8.500377954563628
			BD2G55980	TGACCGAAGGAC	-8.135685675800387
			HighestOrthologScore: -8.135685675800387
			Least Degenerate Matching motifs (multispecies match):
			TGAKTGAYSGAS TGASYGAMSGAC TGAGWGAYSGAS TGAKYGAYSGAS 
			*/
		}
		
		
		for (Map.Entry<String,Set<String>> e : classGeneMap.entrySet()){
			System.out.println(e.getKey()+"\t"+e.getValue().size());
			
		}
		System.out.println("");
		System.out.println("");
		for (Map.Entry<String,Set<String>> e : classGeneMap.entrySet()){
			System.out.println(e.getKey()+"\t"+e.getValue().size());
			
			for (String s : e.getValue()){
				System.out.println(s);
			}
			System.out.println("");
			System.out.println("");
		}
		
	}
	
	private static Scanner createScanner(String line){
		Scanner s = new Scanner(line);
		s.useLocale(Locale.US);
		return s;
	}
	

	private static Map<String, Set<String>> initializeClassGeneMap() {
		
		//classification:
		//Bolduc?				By - Bn
		//OrthologMatch?		Oy - On
		//PWMSC maize >=0		SMy - SMn
		//PMWSC Orthologs >=0   SOy - SOn
		
		
		Map<String, Set<String>> m = new HashMap<String, Set<String>>();
		
		for (int i=0; i<2; i++){
			for (int j=0; j<2; j++){
				for (int k=0; k<2; k++){
					for (int l=0; l<2; l++){
						m.put(generateKeyForClass(i==0,j==0,k==0,l==0), new HashSet<String>());
					}
				}
			}
		}
		
		return m;
		
		
	}
	
	private static String generateKeyForClass(boolean b, boolean o, boolean pm, boolean po){
		char yes = 'y';
		char no = 'n';
		
		return  (Bolduc + (b?yes:no) + OrthoMatch + (o?yes:no) 
					+ PWMSCMaize + (pm?yes:no) + PWMSCOrtho +(po?yes:no)); 
		
		
		
	}
	
}
