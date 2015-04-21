package be.iminds.cloudspeller.PatternMatching;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;
import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.MotifFactory;
import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.output.MotifBLSRestrictions;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.LineIterator;

public class OCRPreparator {

	static {
		int [] thresholds = {15,50,60,70,90,95};
		BLS.initializeBLSConstants(thresholds);
		FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
	}
	
	
	public static void main(String [] args) throws IOException{
		
		MotifBLSRestrictions restrictions = new MotifBLSRestrictions(0,50);
		
		//String exp = "C99F500BLS50";
		String confChartsFile = "/home/ddewitte/Desktop/CloudSpellerExperimentsFinal/PermGroupBestMotifs/BLS95/charsKBestBLS95.txt";
		                                                                             
		//String confChartsFile = "F500C90BLS950All.txt";

		String pmFile = "/home/ddewitte/Desktop/PatternsKBestBLS95.txt";
	
		generatePMFile(confChartsFile,pmFile,restrictions);
		
	}

	private static void generatePMFile(String confChartsFile, String pmFile,
			MotifBLSRestrictions restrictions) throws IOException {
		
		LineIterator iterator = new LineIterator(confChartsFile);
		MotifFactory factory = new IUPACFactory(IUPACType.FULL);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(pmFile));
		
		while (iterator.hasNext()){
			String line = iterator.next();
			
			BLSConfidenceGraph graph = new BLSConfidenceGraph(line,factory);
			
			boolean [] restr = restrictions.checkRestrictions(graph);
			
			//System.out.println(graph);
			
			for (int i=0; i<restr.length; i++){
				if (restr[i]){
					//System.out.println(graph.getMotif());
					writer.append(graph.getMotif().toString());
					writer.append('\t');
					writer.append(""+BLS.getBLSThresholds()[i]);
					writer.newLine();
					break;
				}
			}
			writer.flush();
			
		}
		
		writer.close();
		
		
	}
	
}
