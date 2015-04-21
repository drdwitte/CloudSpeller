package be.iminds.cloudspeller.input;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import be.iminds.cloudspeller.toolbox.LineIterator;

public class ShuffledABDatasetPreparator {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		String outputDir="/home/ddewitte/Desktop/FilesForABonVW/";
		
		String filename = "/home/ddewitte/Desktop/groupOrthoAll.txt";
		
		LineIterator iterator = new LineIterator(filename);
		
		BufferedWriter writer = null;
		
		while (iterator.hasNext()){
			
			String line = iterator.next();
			if (line.length()==0){
				continue;
			}
			
			if (line.startsWith("iORTHO")){
				if (writer!=null)
						writer.close();
				
				String outputFile = line.trim()+".txt";
				writer = new BufferedWriter(new FileWriter(outputDir+outputFile));
				
				/*String temp=*/ iterator.next();
				int numberOfGenes = Integer.parseInt(iterator.next());
				
				for (int i=0; i<numberOfGenes; i++){
					line = iterator.next();
					writer.append(">"+line.split("\t")[0]);
					writer.newLine();
					writer.append(iterator.next());
					writer.newLine();
					writer.newLine();
				}
			}
		}
	
		writer.close();
		
		
		


		
	}
	
	
	

}
