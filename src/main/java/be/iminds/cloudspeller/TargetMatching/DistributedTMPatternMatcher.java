package be.iminds.cloudspeller.TargetMatching;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.util.ToolRunner;

import be.iminds.cloudspeller.driver.DistributedPatternMatcher;


public class DistributedTMPatternMatcher extends DistributedPatternMatcher {
	
	@Override
	public int run(String[] args) throws Exception {
 		
		//-patterns patternsDIR be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/
		if (args.length!=4){
			printUsage();
			return 1;
		}
		
		List<String> files = getFiles(args[1]);
		String inputDataset=args[2];
		String outputDir=args[3];
		if (outputDir.endsWith("/")){
			
		} else {
			outputDir+="/";
		}
		
		
		for (int i=0; i<files.size(); i++){
			
			
			String [] newArgs= new String[4];
			newArgs[0]="-patterns";
			newArgs[1]=files.get(i);
			newArgs[2]=inputDataset;
			
			String tmFile = files.get(i);
			int start=tmFile.lastIndexOf("/")+1;
			int stop = tmFile.lastIndexOf(".");
			String motifName = tmFile.substring(start,stop);
			
			newArgs[3]=outputDir+motifName+"/";
			System.err.println("be.iminds.cloudspeller.output: "+newArgs[3]);
			
			
			super.run(newArgs);
		}

		return 0;
	}
	
	private List<String> getFiles(String patternsDIR) throws IOException {
		List<String> filenames = new ArrayList<String>();
		FileSystem fs = FileSystem.get(getConf());
        FileStatus[] status = fs.listStatus(new Path(patternsDIR));
        for (int i=0;i<status.length;i++){
        	filenames.add(status[i].getPath().toString());
        }
        return filenames;
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DistributedTMPatternMatcher(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
		
		
	protected static int printUsage() {
		System.out.println("bin/hadoop jar DBTM.jar -patterns patternsDIR be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
	}	

}
