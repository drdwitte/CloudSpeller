package be.iminds.cloudspeller.postprocessing_Single;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;


import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;



public class AdaptorPvalueToScorer {

	static final int [] t = {15,50,60,70,90,95};
	
	public static void main (String [] args) throws IOException{
		
		BLS.initializeBLSConstants(t);
		FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
		String dir = "/home/ddewitte/Bureaublad/TestMotifHadoop/outputPvalueSim/";
		String filename = "part-r-00000";
		String output = "PermDistr.txt";
		
		BufferedReader in = new BufferedReader(new FileReader(dir+filename));
		BufferedWriter out =new BufferedWriter(new FileWriter(dir+output));
		String line;
		
		while ( (line=in.readLine()) !=null){
			if (line.length()==0)
				continue;
			
			String permKey=line.trim();
			line=in.readLine(); //skip
			
			
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				line = in.readLine();
				Scanner scanner = GeneralToolbox.generateScanner(line);
								
				int bls = scanner.nextInt();
				double mu = scanner.nextDouble();
				double sigma = scanner.nextDouble();
				
				out.write(permKey+"_"+bls+"\t"+mu+"\t"+sigma+"\n");
				
				scanner.close();
			}
			
			line = in.readLine(); //skip
			
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				line = in.readLine(); //skip normality tests
			}
			
			
			
			
		}
		
		
		in.close();
		
		
	}
}
