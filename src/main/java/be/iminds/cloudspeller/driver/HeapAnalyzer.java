package be.iminds.cloudspeller.driver;

public class HeapAnalyzer {

	private static final int mb = 1024*1024;
	private Runtime runtime;
	
	
	public HeapAnalyzer(){
		runtime = Runtime.getRuntime();
	}
	
	public void printUsedMemory(){
		System.err.println("Used Memory:"
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
	}
	
	public void printAvailableMemory(){
		System.err.println("Available Memory:"
	            + runtime.totalMemory() / mb);
	}
	
	public void printMaximalAvailableMemory(){
		System.err.println("Max Memory:" + runtime.maxMemory() / mb);
	}
	
	public void printAll(String message){
		System.err.println(message);
		printAllNoMessage();
	}
	
	public void printAllNoMessage(){
		printUsedMemory();
		printAvailableMemory();
	}
}
