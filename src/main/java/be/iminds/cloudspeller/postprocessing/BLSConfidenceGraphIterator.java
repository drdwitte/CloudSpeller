package be.iminds.cloudspeller.postprocessing;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;

public class BLSConfidenceGraphIterator implements Iterator<BLSConfidenceGraph> {

	private BufferedReader in;
	private String line;
	private MotifFactory factory = new IUPACFactory(IUPACType.FULL);
	private boolean terminated=false;
	
	
	public BLSConfidenceGraphIterator(String filename) throws IOException{
	
		in = new BufferedReader(new FileReader(filename));
		line = in.readLine();
	}
	
	@Override
	public boolean hasNext()  {
		boolean hn = line!=null; 
		if (!hn){
			terminate();
		}
		return hn;
	}

	@Override
	public BLSConfidenceGraph next() {
		BLSConfidenceGraph graph = new BLSConfidenceGraph(line,factory);
		try {
			line = in.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return graph;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	
	
	public void terminate(){
		if (terminated){ return; } //trying to close it twice!
				
		try {
			in.close();
			terminated=true;
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
