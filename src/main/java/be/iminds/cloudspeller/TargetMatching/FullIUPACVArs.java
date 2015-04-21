package be.iminds.cloudspeller.TargetMatching;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import be.iminds.cloudspeller.KN1Analysis.DegMatchesProcessor;
import be.iminds.cloudspeller.KN1Analysis.KN1Toolbox;
import be.iminds.cloudspeller.KN1Analysis.ScorePair;

public class FullIUPACVArs {

	public static final double minimumPWMScore=11;
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		KN1Toolbox.setKN1Directory("/home/ddewitte/Bureaublad/CloudSpellerExperimentsFinal/KN1_6NOV/");
		String dir = KN1Toolbox.getKN1Directory();
		
		String regex = "TGA..GA..GA.";
		List<StringBuilder> matches = DegMatchesProcessor.createAllDegMatchesWithRegExp(regex, "ACGTNMRSWYKBDHV");
		System.out.println("num matches: "+matches.size());
		
		
		DegMatchesProcessor processor = new DegMatchesProcessor(null);
		//processor.initializeBSMap(regex,KN1Toolbox.generateModifiedPWM());
		processor.initializeBSMap(regex,KN1Toolbox.generateKN1PWM());

		
		String patternOutput = "PatternsFullIUPAC.txt"; //PWMMOD.txt";
		String recordOutput = "RecordsFullIUPAC.txt";   //PWMMOD.txt";
		
		BufferedWriter outp = new BufferedWriter(new FileWriter(dir+patternOutput));
		BufferedWriter outr = new BufferedWriter(new FileWriter(dir+recordOutput));
		
		int numberOfFullIUPACOnlies=0;
		int numberWithPosScore=0;
		
		int blsMin = 15;
		
		
		SortedSet<ScorePair> scoreTableWorstBS = new TreeSet<ScorePair>();

		
		for (StringBuilder m : matches){
			
			int i1 = m.indexOf("B");
			int i2 = m.indexOf("D");
			int i3 = m.indexOf("H");
			int i4 = m.indexOf("V");
			
			int sum = i1+i2+i3+i4;
			
			if (sum==-4){
				continue;
			} else {
				//System.out.println(m);
				String motif = m.toString();
				Double score = processor.calculatescoreWorstBS(motif);
				numberOfFullIUPACOnlies++;
				if (score >= minimumPWMScore){
				
					numberWithPosScore++;
					scoreTableWorstBS.add(new ScorePair(motif,score));
					outp.write(motif+"\t"+blsMin+"\n");
				}
				
				
			}
			
			
		}
		
		outp.close();
		
		
		for (ScorePair p : scoreTableWorstBS){
		
			outr.write(p.getKey()+"\t15\t3\t100.0\t"+p.getValue()+"\n");
		}
		outr.close();
		
		System.out.println("Number of iupac onlies "+numberOfFullIUPACOnlies);
		System.out.println("Number of iupac onlies filtered: "+numberWithPosScore);
		
	}

}
