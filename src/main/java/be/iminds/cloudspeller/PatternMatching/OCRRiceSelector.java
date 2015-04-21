package be.iminds.cloudspeller.PatternMatching;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import be.iminds.cloudspeller.toolbox.LineIterator;

public class OCRRiceSelector {

	public static void main (String [] args) throws IOException{
		
		
		String file = args[0];
		
		String output = file+"Filtered"+".txt";
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		
		LineIterator iterator = new LineIterator(file);
		
		pmlines:while (iterator.hasNext()){
			String line = iterator.next();
			String [] splits = line.split("\t");
			
			
			
			
			if (splits[2].startsWith("OS")){
				
				if (splits[0].startsWith("N")){ 
					continue pmlines;
				}
				
				if (splits[0].split("_")[0].endsWith("N")){
					continue pmlines;
				}
								
				writer.append(line);
				writer.newLine();
				writer.flush();
			}
		}
		
		writer.close();
	}
	
	
}
