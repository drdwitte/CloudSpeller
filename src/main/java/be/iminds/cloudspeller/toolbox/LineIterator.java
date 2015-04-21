package be.iminds.cloudspeller.toolbox;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;


public class LineIterator implements Iterator<String> {

	private BufferedReader in;
	private String line;
	private boolean terminated=false;
	
	public LineIterator(String filename) throws IOException {
		
		in = new BufferedReader(new FileReader(filename));
		line = in.readLine();
	}
	
	@Override
	public boolean hasNext() {
		boolean hn = (line!=null); 
		
		if (!hn){
			terminate();
		}
		return hn;
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
	@Override
	public String next() {
		String oldLine = new String(line);		
		try {
			line = in.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return oldLine;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	} 

}
