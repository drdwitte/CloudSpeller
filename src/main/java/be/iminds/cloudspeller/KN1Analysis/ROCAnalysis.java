package be.iminds.cloudspeller.KN1Analysis;


import java.io.IOException;



public class ROCAnalysis {

	
	public static void main(String [] args) throws IOException{
		//filter deNovoRecords
		
		/*String exp = "C90";
		String inputDir = "/home/ddewitte/Bureaublad/CloudSpellerExperimentsFinal/KN1_6NOV/DeNovoMatchesFiltered/";
		String inputFile = inputDir+"deNovoMatchesFiltered"+exp+"F0.txt";
		int FThreshold=50;
		
		
		
		String outputFile = inputDir+"deNovoMatchesFiltered"+exp+"F"+FThreshold+".txt";
		LineIterator iterator = new LineIterator(inputFile);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		
		while (iterator.hasNext()){
			
			
			DeNovoRecord record = new DeNovoRecord(iterator.next());
			
			int F = record.getNumFamilies();
			
			if (F>=FThreshold){
				writer.write(record.toString()); writer.newLine();
			}
		}
		
		writer.close();
		*/
		
		KN1Toolbox.setExperiment("C90F0");
		String outputPM = KN1Toolbox.getPMOutput(KN1Toolbox.getExperiment());
				
		KN1Toolbox.setExperiment("C90F50");
		String recordsFile = KN1Toolbox.getDeNovoRecords(KN1Toolbox.getExperiment());
		DegMatchesProcessor.reworkPMOutput(recordsFile, outputPM,KN1Toolbox.getExperiment());

		
	}
}
