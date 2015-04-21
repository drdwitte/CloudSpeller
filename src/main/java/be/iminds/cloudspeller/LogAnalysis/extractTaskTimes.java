package be.iminds.cloudspeller.LogAnalysis;

import java.io.IOException;

import be.iminds.cloudspeller.toolbox.LineIterator;

public class extractTaskTimes {
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		String dir = "/home/ddewitte/Bureaublad/Back2Amazon/BigSimulationInfo/";
		String filename = "rawCopyPasteData.txt"; 
				
		
		
		LineIterator iterator = new LineIterator(dir+filename);
		
		while (iterator.hasNext()){
			String line = iterator.next();
			
			int start = line.indexOf('(');
			int stop = line.indexOf(')');
			//String mapTime;
			if (start >= 0 && stop >= 0){
				String maptime = line.substring(start+1,stop);
				
				String [] splits = maptime.split(",");
				//System.out.println(maptime);
				
				int hours=0; int mins=0; int secs=0;
				
				
				for (int i=0; i<splits.length; i++){
				
					//System.out.println(splits[i]);
					String sub;
					if (splits[i].contains("h")){
						sub=splits[i].substring(0,splits[i].indexOf("h")).trim();
						if (sub.length()>0){
							hours = Integer.parseInt(sub);
						}
					} else if (splits[i].contains("m")){
						sub=splits[i].substring(0,splits[i].indexOf("m")).trim();
						if (sub.length()>0){
							mins = Integer.parseInt(sub);
						}
					} else if (splits[i].contains("s")){
						sub=splits[i].substring(0,splits[i].indexOf("s")).trim();
						if (sub.length()>0){
							secs = Integer.parseInt(sub);
						}
					}
					
				}

				//System.out.println(hours+":"+mins+":"+secs);
				System.out.println(3600*hours+60*mins+secs);
			}
			
			
		}
		
		
		
		
		
		
		
		
		
		
		
	}

}
