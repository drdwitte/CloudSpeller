package be.iminds.cloudspeller.output;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.Motif;

public class BLSConfGraphFactory implements ConfidenceGraphFactory{

	@Override
	public ConfidenceGraph createGraph(Motif motif, FreqVec freqVec,
			ProbabilityVector probVec) {
		return new BLSConfidenceGraph(motif, freqVec, probVec);
	}
	

	/**
	 * Read be.iminds.cloudspeller.output of the following form:
	 */
	/*
	CCCT	Motif=CCCT
	40	0.0	44	
	50	0.0	44	
	60	0.0	44	
	70	0.0	44	
	80	0.0	40	
	90	0.0	38
	*/
	public void getGraphsFromBuffer(BufferedReader in, Map<Motif,BLSConfidenceGraph> graphs) throws IOException{
		
		
		
		IUPACFactory fac = new IUPACFactory(IUPACType.FULL);
		String line;
		while ((line=in.readLine())!=null){
			
			Motif motif;
			FreqVec freqVec;
			ProbabilityVector probVec;
			
			if (line.length()==0){
				continue; //skip empyty lines (null is only of EOF)
			}
			
			Scanner scan=new Scanner(line);
			
			scan.next(); //permutation group ID
			String motifLine = scan.next();
			scan.close();
			String m = motifLine.split("=")[1];
			
			
		
			motif = fac.createMotifFromString(m);
			freqVec = new FreqVec();
			probVec = new ProbabilityVector();
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				line = in.readLine();
				scan = new Scanner(line);
				scan.useLocale(Locale.US); //interprete . as decimal sign -> standard is , !!
				scan.nextInt(); //ignore BLS
				probVec.setProb(i,scan.nextDouble());
				freqVec.setFreq(i,scan.nextInt());
				scan.close();
			}
			
			BLSConfidenceGraph output = new BLSConfidenceGraph(motif, freqVec, probVec);
			graphs.put(motif,output);
			
			if (graphs.size()%10000==0){
				System.out.println(graphs.size()+" aws graphs read");
			}
			
		}
	}


	public void getFreqVecsFromBuffer(BufferedReader in,
			Map<Motif, FreqVec> awsOutput) throws IOException {
		
		
		IUPACFactory fac = new IUPACFactory(IUPACType.FULL);
		String line;
		while ((line=in.readLine())!=null){
		
			if (line.length()==0){
				continue; //skip empyty lines (null is only of EOF)
			}
			
			Scanner scan=new Scanner(line);
			
			scan.next(); //permutation group ID
			String motifLine = scan.next();
			scan.close();
			String m = motifLine.split("=")[1];
			
			Motif motif = fac.createMotifFromString(m);
			
			FreqVec freqVec = new FreqVec();
			
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				line = in.readLine();
				scan = new Scanner(line);
				scan.useLocale(Locale.US); //interprete . as decimal sign -> standard is , !!
				scan.nextInt(); //ignore BLS
				scan.nextDouble();
				freqVec.setFreq(i,scan.nextInt());
				scan.close();
			}
			
			
			awsOutput.put(motif,freqVec);
			
			
			if (awsOutput.size()%10000==0){
				System.out.println(awsOutput.size()+" aws graphs read");
			}
			
		}
	}

}
