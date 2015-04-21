package be.iminds.cloudspeller.postprocessing_Single;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import be.iminds.cloudspeller.driver.MotifReferenceDistanceCalculator;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.postprocessing.CharacterDistanceCalculator;
import be.iminds.cloudspeller.postprocessing.MotifAligner;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;


public class DBWangMatcher extends Configured implements Tool {

	protected String input;
	protected String output;
	protected static final String paste="_";
	protected static final double smallDouble = 1e-6;
	
	/**
	 * No motif filtering! Use DBMotifFilter for this!
	 * @author ddewitte
	 *
	 */
	public static class WangMapper extends Mapper<LongWritable, Text, Text,Text> {
		
		private static Set<String> wangMotifs;
		private static Text outputKey = new Text();
		private static Text outputValue = new Text();
		//private static ConfidenceGraphRestrictions restrictions;
		private static Map<String,Double> wangDistance = new HashMap<String,Double>();
		private static Map<String,Set<String>> wangMatches = new HashMap<String,Set<String>>();
		//private static final int [] t = {15,50,60,70,90,95};
		//private static MotifFactory motifFactory;

		private static Alphabet alphabet;
		private static CharacterDistanceCalculator calculator; 
		private static MotifAligner aligner;
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			
			alphabet = new IUPACAlphabet(IUPACType.TWOFOLDSANDN);
			
			calculator = new MotifReferenceDistanceCalculator(alphabet);  
			aligner = new MotifAligner(alphabet, calculator);

			wangMotifs = new HashSet<String>();
			
			Path[] uris = DistributedCache.getLocalCacheFiles(context
				    .getConfiguration());
			
			String wangFile = uris[0].toString();
	    
			
		    BufferedReader in = new BufferedReader(new FileReader(wangFile));
		    String wangMotifLine;
		    while ((wangMotifLine = in.readLine()) != null) {
		    	
		    	wangMotifs.add(wangMotifLine.split("\t")[0]);
		    }
		    
		    in.close();
		    
		    //initialize wang map
			
			
			for (String s : wangMotifs){
				wangDistance.put(s,Double.MAX_VALUE);
				wangMatches.put(s,new HashSet<String>());
			}
		    
			
		}
		
		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			
			for (String wang : wangMotifs){
				
				double d = wangDistance.get(wang);
				outputKey.set(wang);
								
				for (String match : wangMatches.get(wang)){
					outputValue.set(match+paste+d);
					context.write(outputKey,outputValue);
				}
				
				
			}
			
			wangMotifs.clear();
			wangMatches.clear();
			wangDistance.clear();
			wangMotifs=null;
			wangMatches=null;
			wangDistance=null;
			
		}
		
		
		@Override
		protected void map(LongWritable key, Text value, Context context)
		throws IOException, InterruptedException {
			
			String line = value.toString();
			if (line.length()==0)
				return;
			
			Scanner scanner = GeneralToolbox.generateScanner(line);
			scanner.next(); //key
			String m = scanner.next(); //motif

			IUPACMotif motif = new IUPACMotif(m,0);
			String  mRev = motif.getComplement().toString();
		
			
			System.out.println(m+"\t"+mRev);
			
			for (String wang : wangMotifs){
			
			
				double d1=calculateDistance(wang,m);
				double d2=calculateDistance(wang,mRev);
				double dNew = Math.min(d1,d2);
				
				double dCurrent = wangDistance.get(wang);
				
				if (smallerThan(dNew,dCurrent)){
					wangDistance.put(wang,dNew);
					wangMatches.get(wang).clear();
					wangMatches.get(wang).add(m);
					wangMatches.get(wang).add(mRev);
				} else if (almostEqual(dNew,dCurrent)){
					wangMatches.get(wang).add(m);
					wangMatches.get(wang).add(mRev);
				} else {
					//do nothing
				}
			
			}
			
			//emit best results
			
			
			
			for (String wang : wangMotifs){
				
				outputKey.set(wang);
				
				double dmin = wangDistance.get(wang);
				
				for (String match : wangMatches.get(wang)){
					outputValue.set(match+paste+dmin);
					context.write(outputKey,outputValue);
				}
				
				
				//outputValue.set();
				
				
			}
			
			
			
		}
		
		double calculateDistance(String wangMotif, String foundMotif){
			return aligner.calculateSimilarityDistance(wangMotif,foundMotif);
		}
			
			
	
				
				
			
	}
	
	public static boolean almostEqual(double d1, double d2){
		
		return Math.abs(d1-d2)<smallDouble;
		
	}
	
	public static boolean smallerThan(double d1, double d2){
		if (!almostEqual(d1,d2)){
			return d1<d2;
			
		} else 
			return false;
	}
	
	
	public static class WangReducer extends Reducer<Text, Text, Text,Text>{
		
		private static Text outputKey;
		private static Text outputValue = new Text();
		
		@Override
		protected void reduce(Text wangMotif, Iterable<Text> motifDistancePairs,
				Context context)
				throws IOException, InterruptedException {
			
			
			outputKey = wangMotif;
		
			double wangDistanceMininum = Double.MAX_VALUE;
			Set<String> wangMatches = new HashSet<String>();
			
			for (Text pair : motifDistancePairs){
				String [] spl = pair.toString().split(paste);
				String motif = spl[0];
				double d = Double.parseDouble(spl[1]);
				
				if (smallerThan(d,wangDistanceMininum)){
					wangDistanceMininum=d;
					wangMatches.clear();
					wangMatches.add(motif);
				} else if (almostEqual(d,wangDistanceMininum)){
					wangMatches.add(motif);
				} else {
					//do nothing
				}
			}
			
			for (String match : wangMatches){
				
				outputValue.set(match+paste+wangDistanceMininum);
				context.write(outputKey,outputValue);
			}
			
		}

		
		
	}
	
	
	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DBWangMatcher(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
	
	
	@Override
	public int run(String[] args) throws Exception {
        
        //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "WangMatcher");
		
		job.setJarByClass(DBWangMatcher.class);
		
  	  	//SETTING MAPREDUCE CLASSES
 	  	
		// IO
		job.setInputFormatClass(TextInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(Text.class);
		job.setMapOutputValueClass(Text.class);

		// (K3,V3)
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(Text.class);

		// set hadoop methods
		job.setMapperClass(WangMapper.class);
		job.setReducerClass(WangReducer.class); 

		processCommandLineArgs(args, job);
		
		try {
			FileInputFormat.setInputPaths(job, new Path(input));
			FileOutputFormat.setOutputPath(job, new Path(output));
		} catch (Exception e) {
			System.err.println(e.getMessage());
		}

		return job.waitForCompletion(true) ? 0 : 1;
	}
	
	protected void processCommandLineArgs(String[] args, Job job) {

		// process command line arguments
		List<String> other_args = new ArrayList<String>();
		
		String patternsFile = "";

		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-r".equals(args[i])) {
					job.setNumReduceTasks(Integer.parseInt(args[++i]));
				} else if ("-wang".equals(args[i])){
					patternsFile=args[++i];
				} else {
					other_args.add(args[i]);
				}
			} catch (NumberFormatException except) {
				System.out.println("ERROR: Integer expected instead of " + args[i]);
				printUsage();
			} catch (ArrayIndexOutOfBoundsException except) {
				System.out.println("ERROR: Required parameter missing from " + args[i - 1]);
				printUsage();
			}
			
		}

		// Make sure there are exactly 2 parameters left.
		if (other_args.size() != 2) {
			System.out.println("ERROR: Wrong number of parameters: " + other_args.size() + " instead of 2.");
			printUsage();
		}
		
		//settingsFile = other_args.get(0);
		input = other_args.get(0);
		output = other_args.get(1);
		
		try {
			DistributedCache.addCacheFile(new URI(patternsFile), job.getConfiguration());
		
		} catch (URISyntaxException e) {
			System.out.println("URI exception: "+patternsFile);
			e.printStackTrace();
		}
	    

	}


	protected static int printUsage() {
		System.out.println("bin/hadoop jar WangMatcher.jar -wang wang.txt be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
	}

	
	
}
