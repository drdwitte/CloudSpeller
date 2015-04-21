package be.iminds.cloudspeller.postprocessing_Single;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;

import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import org.apache.hadoop.mapreduce.Mapper;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.SimultaneousOccurrenceAndConfidenceFiltering;
import be.iminds.cloudspeller.phylogenetics.BLS;

/**
 * Motif -> (prefix,Motif) -> (prefix, #unique motifs) 
 * (important since sometimes motifs and reverse motifs are in be.iminds.cloudspeller.output)
 */


public class DBMotifCounter extends Configured implements Tool {

	private static String input;
	private static String output;
	
	public static class MotifMapper extends Mapper<LongWritable,Text,Text,Text> {

		private static MotifFactory motifFactory;
		private static final int [] t = {15,50,60,70,90,95};
		
		private static Text outputKey = new Text();
		private static Text outputValue = new Text();
		private static final int prefixLength = 5;
		
		private static ConfidenceGraphRestrictions restrictions;
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			restrictions = new SimultaneousOccurrenceAndConfidenceFiltering(localConfCutoff,localFamCutoff);
		}

		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			String strValue = value.toString();
			if (strValue.length()==0)
				return;
			
			BLSConfidenceGraph graph = new BLSConfidenceGraph(strValue, motifFactory);
		
			if (!restrictions.checkRestrictions(graph)){
				return;
			}
			
			Motif m = graph.getMotif();
			Motif mRev = m.getComplement();
			
			//System.out.println(m + "\t -> \t"+mRev);
			
			outputKey.set(m.toString().substring(0,prefixLength));
			outputValue.set(m.toString());
			context.write(outputKey,outputValue);
						
			outputKey.set(mRev.toString().substring(0,prefixLength));
			outputValue.set(mRev.toString());
			context.write(outputKey,outputValue);
		}
		
	}
	
	static class MotifReducer extends Reducer<Text,Text,Text,IntWritable> {

		private static Text outputKey = new Text();
		private static IntWritable outputValue = new IntWritable();
		
		private static Set<String> uniqueMotifs;
		
		@Override
		protected void cleanup(Context context)
				throws IOException, InterruptedException {
			uniqueMotifs = null;
		}


		@Override
		protected void reduce(Text key, Iterable<Text> values,
				Context context)
				throws IOException, InterruptedException {
			
			uniqueMotifs = new HashSet<String>();
			outputKey  = key;
			for (Text t : values){
				uniqueMotifs.add(t.toString());
			}
			outputValue.set(uniqueMotifs.size());
			
			context.write(outputKey,outputValue);
			uniqueMotifs.clear(); uniqueMotifs=null;
		}



		
	}
	
	
	@Override
	public int run(String[] args) throws Exception {
		
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "MotifCounter");
		
		job.setJarByClass(DBMotifCounter.class);
		
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
		job.setMapperClass(MotifMapper.class);
		job.setReducerClass(MotifReducer.class);

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
		System.out.println("bin/hadoop jar MotifCounter -C 50 -F 3 be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DBMotifCounter(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
	
	

}