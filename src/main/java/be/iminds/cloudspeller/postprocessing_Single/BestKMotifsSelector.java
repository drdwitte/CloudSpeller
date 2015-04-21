package be.iminds.cloudspeller.postprocessing_Single;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
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
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.MotifBLSRestrictionsWithCutoff;
import be.iminds.cloudspeller.output.SimultaneousOccurrenceConfidenceAndBLSFiltering;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;

public class BestKMotifsSelector extends Configured implements Tool {

	private static String input;
	private static String output;

	
	public static class MotifPermGroupMapperInMapperComb extends Mapper<LongWritable,Text,Text,Text> {

		private static final int [] t = {15,50,60,70,90,95};

		private static int counter=0;
		private static long start;
 		private static int k;

		private static MotifBLSRestrictionsWithCutoff restrictions;
		private static Text mapOutputKey = new Text();
		private static Text mapOutputValue = new Text();
		
		private static String currentKey = "";
		private static List<String> graphs = new ArrayList<String>(50000);
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			start=System.nanoTime();
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			int localBLSCutoff = context.getConfiguration().getInt("BLS",0);
			
			k = context.getConfiguration().getInt("KBest",10);

			restrictions = new MotifBLSRestrictionsWithCutoff(localConfCutoff, localFamCutoff, localBLSCutoff );
		}
		
		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			++counter;
				
			
			String strValue = value.toString();
			if (strValue.length()==0)
				return;
			
			restrictions.getMaxFamilies(strValue);
	
			Scanner scanner = GeneralToolbox.generateScanner(strValue);
		
			String localKey=scanner.next();
			
			
			
			if (localKey.equals(currentKey)){
				graphs.add(strValue);
				
			} else {
			
				
				System.out.println("S before flush: "+graphs.size()+"\t tot: "+counter + "\t"+(System.nanoTime()-start)/Math.pow(10,9));
			
				flushGraphs(context);
				currentKey=localKey;
				graphs.add(strValue);
			}
		
			

		}
		
		private void flushGraphs(Context context) throws IOException, InterruptedException {
			
			//determine top k graphs
			
			BestMotifContainer bestMotifs = new BestMotifContainer(k,restrictions);
			for (String g : graphs){
				
				bestMotifs.add(g);
				
			}

			List<String> l = bestMotifs.getAll();
			mapOutputKey.set(currentKey);
			for (String g : l){
				mapOutputValue.set(g);
				
				context.write(mapOutputKey,mapOutputValue);
				
				System.out.println("write: "+mapOutputKey+"\t"+mapOutputValue);
			}
			
			


			
			
			graphs.clear();

			
		}

		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			flushGraphs(context);
		}
		
		
	}
	
	
	public static class MotifPermGroupMapper extends Mapper<LongWritable,Text,Text,Text> {

		private static MotifFactory motifFactory;
		private static final int [] t = {15,50,60,70,90,95};
		
				
		private static ConfidenceGraphRestrictions restrictions;
		private static Text mapOutputKey = new Text();
		private static Text mapOutputValue = new Text();
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			int localBLSCutoff = context.getConfiguration().getInt("BLS",0);
			
			restrictions = new SimultaneousOccurrenceConfidenceAndBLSFiltering(localConfCutoff, localFamCutoff, localBLSCutoff );
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
			
			Scanner s = GeneralToolbox.generateScanner(strValue);
			mapOutputKey.set(s.next());
			s.close();
			
			graph.toText(mapOutputValue);
			context.write(mapOutputKey,mapOutputValue);

		}
		
	}
	
 	public static class MotifPermGroupReducer extends Reducer<Text,Text, NullWritable, Text> {

		private static final int [] t = {15,50,60,70,90,95};

		private static MotifBLSRestrictionsWithCutoff restrictions;
		private static Text outputValue = new Text();
			
		
 		private static int k;
 	
 		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
 			
 			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			int localBLSCutoff = context.getConfiguration().getInt("BLS",0);
			k = context.getConfiguration().getInt("KBest",10);
 			
			restrictions = new MotifBLSRestrictionsWithCutoff(localConfCutoff, localFamCutoff, localBLSCutoff );
 		}


	
		@Override
		protected void reduce(Text key, Iterable<Text> values, Context context)
				throws IOException, InterruptedException {
			
			BestMotifContainer bestMotifs = new BestMotifContainer(k,restrictions);
			
			for (Text value : values) {
							
				bestMotifs.add(value.toString());
			}
	
			List<String> l = bestMotifs.getAll();
			
			for (String g : l){
				outputValue.set(g);
				context.write(NullWritable.get(),outputValue);
			}
	
			
		}
 	}

 	
 	
	@Override
	public int run(String[] args) throws Exception {
		
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "BestKMotifsSelector");
		
		job.setJarByClass(BestKMotifsSelector.class);
		
   	  	//SETTING MAPREDUCE CLASSES
  	  	
		// IO
		job.setInputFormatClass(TextInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(Text.class);
		job.setMapOutputValueClass(Text.class);

		//(K3,V3)
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(Text.class);

		// set hadoop methods
		job.setMapperClass(MotifPermGroupMapperInMapperComb.class);
		job.setReducerClass(MotifPermGroupReducer.class);
		
		
		processCommandLineArgs(args, job);
		
		// for multiple runs remove be.iminds.cloudspeller.output folder of previous run
		//FileSystem fs = FileSystem.get(job.getConfiguration());
		//fs.delete(new Path(be.iminds.cloudspeller.output), true); // true stands for recursively deleting the folder you gave

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
		int minBLS=0;
		int k=10;
		
		List<String> other_args = new ArrayList<String>();
		
		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-C".equals(args[i])) {
					confidenceCutoff = Integer.parseInt(args[++i]);
				} else if ("-F".equals(args[i])){
					familyCutoff = Integer.parseInt(args[++i]);
				
				} else if ("-BLS".equals(args[i])){
					minBLS = Integer.parseInt(args[++i]);	
				} else if ("-KBest".equals(args[i])){
					k = Integer.parseInt(args[++i]);		
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
		job.getConfiguration().setInt("BLS",minBLS);
		job.getConfiguration().setInt("KBest",k);
		
		
		if (other_args.size()!=2){
			System.out.println("ERROR: Wrong number of parameters: " + args.length + " instead of 8.");
			printUsage();
		}
		
		
		input = other_args.get(0);
		output = other_args.get(1);
	}

	private static int printUsage() {
		System.out.println("bin/hadoop jar MotifFilter -C 50 -F 3 -BLS -KBest 10 0 be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new BestKMotifsSelector(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
}
