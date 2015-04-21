package be.iminds.cloudspeller.postprocessing;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;

public class ConfGraphCollector {
	
	
	private BLSConfidenceGraph [] graphArray;
	private int currentSize;
	private String output;
	
	
	public ConfGraphCollector(int maxStorage, String output){
		graphArray = new BLSConfidenceGraph[maxStorage];
		currentSize = 0;
		this.output=output;
	}
	
	public void addGraph(BLSConfidenceGraph graph) throws IOException{
		if (currentSize==graphArray.length){
			flushCollector();
		}
		graphArray[currentSize++]=graph;
	}
	
	public void flushCollector() throws IOException{
		boolean append=true;
		BufferedWriter out = new BufferedWriter(new FileWriter(output,append));
		
		for (BLSConfidenceGraph g : graphArray){
			out.write(g.getMotif().createContent().toString());
			out.write("\t");
			out.write(g.toOneLineStringFormat());
		}
		out.close();
		
		currentSize=0;
	}
	
}
