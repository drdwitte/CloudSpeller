package be.iminds.cloudspeller.postprocessing_Single;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
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


import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.postprocessing_Single.DBMotifCounter.MotifMapper;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

public class DBMotifScorer extends Configured implements Tool {

	private static String input;
	private static String output;
	
	
	static class MotifMapperScorer extends Mapper<LongWritable,Text,Text,DoubleWritable> {

		private static MotifFactory motifFactory;
		private static final int [] t = {15,50,60,70,90,95};
		
		private static Text outputKey = new Text();
		private static DoubleWritable outputValue = new DoubleWritable();
		
		private static Map<String,Stats> distributionMap = new HashMap<String,Stats>();
		
		
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
			
			initializeDistributionMap(context);
		}

		private void initializeDistributionMap(
				Context context) throws IOException {
						
			Path[] uris = DistributedCache.getLocalCacheFiles(context
				    .getConfiguration());
			
			String file = uris[0].toString();
	    
			
		    BufferedReader in = new BufferedReader(new FileReader(file));
		    String motifStatsString;
		    
		    //be.iminds.cloudspeller.input: motif_BLS mu sigma
		    while ((motifStatsString = in.readLine()) != null) {
		    	Scanner scanner = GeneralToolbox.generateScanner(motifStatsString);
		    	String group  = scanner.next();
		    	double mu = scanner.nextDouble();
		    	double sigma = scanner.nextDouble();
		    	distributionMap.put(group,new Stats(mu,sigma));
		    	scanner.close();
		    }
		    
		    in.close();
		}

		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			BLSConfidenceGraph graph = new BLSConfidenceGraph(value.toString(), motifFactory);
		

			
			IUPACMotif m = (IUPACMotif) graph.getMotif();
			String content = m.createContent().toString();
			
			
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){
				int bls = BLS.getBLSThresholds()[i];
				
				Stats stats = distributionMap.get(content+"_"+bls);
				double score = 0.0;
				if (stats!=null){
					stats.calcZScore(graph.getFreq(i));
				} else {
					System.err.println(content+"_"+bls+" NOT FOUND!");
				}
				
				if (score >= 1.0){
				
					outputKey.set(m+"_"+bls);
					outputValue.set(score);
					context.write(outputKey, outputValue);
				}
			}
		}
		
		class Stats {
			
			private double mu;
			private double sigma;
			private static final double smallestSigma = 1e-9;

			public double getMu() {
				return mu;
			}


			public double getSigma() {
				return sigma;
			}

			
			
			public Stats (double mu, double sigma){
				this.mu=mu;
				this.sigma=sigma;
				
				if (this.sigma < smallestSigma){
					this.sigma=smallestSigma;
				}
				
				
				
			}
			
			public double calcZScore(int occ){
				return (occ-mu)/sigma;
			}
			
		}
		
	}
	
	


		
	

	
	@Override
	public int run(String[] args) throws Exception {
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "MotifScorer");
		
		job.setJarByClass(DBMotifScorer.class);
		
   	  	//SETTING MAPREDUCE CLASSES
  	  	
		// IO
		job.setInputFormatClass(TextInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(Text.class);
		job.setMapOutputValueClass(DoubleWritable.class);

		// set hadoop methods
		job.setMapperClass(MotifMapper.class);
		job.setReducerClass(Reducer.class);

		processCommandLineArgs(args, job);
		
		// for multiple runs remove be.iminds.cloudspeller.output folder of previous run
		FileSystem fs = FileSystem.get(job.getConfiguration());
		fs.delete(new Path(output), true); // true stands for recursively deleting the folder you gave

		// setting IO

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
		
		String file = "";

		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-r".equals(args[i])) {
					job.setNumReduceTasks(Integer.parseInt(args[++i]));
				} else if ("-stats".equals(args[i])){
					file=args[++i];
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
			DistributedCache.addCacheFile(new URI(file), job.getConfiguration());
		
		} catch (URISyntaxException e) {
			System.out.println("URI exception: "+file);
			e.printStackTrace();
		}
	}

	private static int printUsage() {
		System.out.println("bin/hadoop jar SecondarySorter.jar be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DBMotifScorer(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
}
	