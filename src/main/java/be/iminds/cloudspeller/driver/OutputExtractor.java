package be.iminds.cloudspeller.driver;

import org.apache.hadoop.io.Text;

import be.iminds.cloudspeller.output.ConfidenceGraph;


public interface OutputExtractor {
	
	public void extract(ConfidenceGraph graph);
	public void setKeyValueTextContainers(Text outputKey, Text outputValue);
}

