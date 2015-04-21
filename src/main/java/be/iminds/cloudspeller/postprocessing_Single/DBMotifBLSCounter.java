package be.iminds.cloudspeller.postprocessing_Single;

import be.iminds.cloudspeller.indexing.PatternBLSPair;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import org.apache.commons.lang.NotImplementedException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
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
import be.iminds.cloudspeller.output.MotifBLSRestrictions;

import be.iminds.cloudspeller.phylogenetics.BLS;

public class DBMotifBLSCounter extends Configured implements Tool {

	private static String input;
	private static String output;
	
	
	static class MotifBLSMapper extends Mapper<LongWritable,Text,Text,Text> {

		private static MotifFactory motifFactory;
		private static final int [] t = {15,50,60,70,90,95};
		
		private static Text outputKey = new Text();
		private static Text outputValue = new Text();
		private static final int prefixLength = 4;
		
		private static MotifBLSRestrictions restrictions;
	
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			restrictions = new MotifBLSRestrictions(localFamCutoff,localConfCutoff);
		}
		
		
		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			throw new NotImplementedException();
		}

		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			
			BLSConfidenceGraph graph = new BLSConfidenceGraph(value.toString(),motifFactory);
			boolean [] motifBLSCheck = restrictions.checkRestrictions(graph);
			
			Motif m = graph.getMotif();
			Motif mRev = m.getComplement();
			
			String mPref = m.toString().substring(0,prefixLength);
			String mRevPref = mRev.toString().substring(0,prefixLength);
			
			
			
			for (int i=0; i<motifBLSCheck.length; i++){
				
				if (motifBLSCheck[i]){

					int bls = BLS.getBLSThresholds()[i];
					
					PatternBLSPair p = new PatternBLSPair(m.toString(),bls);
					PatternBLSPair pRev = new PatternBLSPair(mRev.toString(),bls);
					
					outputKey.set(mPref);
					outputValue.set(p.toString());
					context.write(outputKey,outputValue);
					
					outputKey.set(mRevPref);
					outputValue.set(pRev.toString());
					context.write(outputKey,outputValue);
					
				} else {
					//skip
				}
			}
		}

		

		
		
	}
	
	static class MotifBLSReducer extends Reducer<Text,Text,Text,IntWritable> {

		private static Text outputKey = new Text();
		private static IntWritable outputValue = new IntWritable();
		
		private static final int [] t = {15,50,60,80,90,95};
		
		private static Map<Integer,Set<PatternBLSPair>> uniqueMotifBLSPairs;
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			
		}
		
		@Override
		protected void cleanup(Context context)
				throws IOException, InterruptedException {
			uniqueMotifBLSPairs = null;
		}


		@Override
		protected void reduce(Text key, Iterable<Text> values,
				Context context)
				throws IOException, InterruptedException {
			
			
			uniqueMotifBLSPairs = new HashMap<Integer, Set<PatternBLSPair>>();
			
			//initialize
			for (int i=0; i<BLS.getNumberOfIntervals(); i++){
				uniqueMotifBLSPairs.put(BLS.getBLSThresholds()[i],new HashSet<PatternBLSPair>());
			}
			
			
			for (Text t : values){
				
				PatternBLSPair p = new PatternBLSPair(t.toString());
				uniqueMotifBLSPairs.get(p.getBls()).add(p);				
			}
			
			
			
			//write
			for (int i=0; i<BLS.getNumberOfIntervals(); i++){
				
				int bls = BLS.getBLSThresholds()[i];
				int numMotifBLS = uniqueMotifBLSPairs.get(bls).size();
				
				outputKey.set(key.toString()+"_"+bls);
				outputValue.set(numMotifBLS);
				context.write(outputKey,outputValue);
			}
			
			
			
			
			uniqueMotifBLSPairs = null;
			
			
		}



		
	}
	
	
	@Override
	public int run(String[] args) throws Exception {
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "MotifBLSCounter");
		
		job.setJarByClass(DBMotifBLSCounter.class);
		
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
		job.setOutputValueClass(IntWritable.class);

		// set hadoop methods
		job.setMapperClass(MotifBLSMapper.class);
		job.setReducerClass(MotifBLSReducer.class);

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

private static void processCommandLineArgs(String[] args, Job job) {
		

		int confidenceCutoff=50;
		int familyCutoff=3;
		
		List<String> other_args = new ArrayList<String>();
		
		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-C".equals(args[i])) {
					confidenceCutoff = Integer.parseInt(args[++i]);
				} else if ("-F".equals(args[i])){
					familyCutoff = Integer.parseInt(args[++i]);
				
				} else {
					other_args.add(args[i]);
				}
			} catch (NumberFormatException except) {
				System.out.println("ERROR: Integer expected instead of " + args[i]);
				printUsage();
			}
		}
		
		job.getConfiguration().setInt("C",confidenceCutoff);
		job.getConfiguration().setInt("F",familyCutoff);
		
		
		if (other_args.size()!=2){
			System.out.println("ERROR: Wrong number of parameters: " + args.length + " instead of 6.");
			printUsage();
		}
		
		
		
		
		input = other_args.get(0);
		output = other_args.get(1);
	}


	private static int printUsage() {
		System.out.println("bin/hadoop jar SecondarySorter.jar be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DBMotifBLSCounter(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
	
	

}
