package be.iminds.cloudspeller.postprocessing_Comparative;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

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
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.output.MotifBLSRestrictions;
import be.iminds.cloudspeller.phylogenetics.BLS;



public class DBTargetComparator extends Configured implements Tool {

	private static String input;
	private static String output;
	
	static class TargetMapper extends Mapper<LongWritable,Text,Text,Text> {

		private static String filenamePrefix;
		private static MotifFactory motifFactory;
		private static final int [] t = {15,50,60,70,90,95};
		private static MotifBLSRestrictions restrictions;
		
		private static Text outputKey = new Text(); //prefix + BLS
		private static Text outputValue = new Text(); //motif_#targets
		private static final int prefixLength = 4;
		
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			restrictions = new MotifBLSRestrictions(localFamCutoff,localConfCutoff);

			
			FileSplit fileSplit = (FileSplit)context.getInputSplit();
			String tempFilename = fileSplit.getPath().getName();
			filenamePrefix = tempFilename.split(".")[0].substring(0,3);
		}
		
		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			String strValue = value.toString();
			if (strValue.length()==0)
				return;
			
			BLSConfidenceGraph graph = new BLSConfidenceGraph(strValue, motifFactory);
		
			boolean [] motifBLSCheck = restrictions.checkRestrictions(graph);
			
			String motif = graph.getMotif().toString();

			
			char paste = '_';
			for (int i=0; i<motifBLSCheck.length; i++){
				
				if (!motifBLSCheck[i]){
					continue;
				}
			
				//prefix$BLS
				outputKey.set(motif.substring(0,prefixLength+BLS.getBLSThresholds()[i]));
				outputValue.set(motif+paste+filenamePrefix+paste+graph.getFreq(i));
				context.write(outputKey, outputValue);
			
				String revComp = graph.getMotif().getComplement().toString();
				outputKey.set(revComp.substring(0,prefixLength+BLS.getBLSThresholds()[i]));
				outputValue.set(revComp+paste+filenamePrefix+paste+graph.getFreq(i));
				context.write(outputKey, outputValue);

				
			}
		}
		
		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			  
			filenamePrefix = null;

		}
		
	}
	
	
	
	
	
	static class TargetReducer extends Reducer<Text,Text,Text,IntWritable> {
	
		private static Text outputKey = new Text();
		private static IntWritable outputValue = new IntWritable();
		private static final String paste = "_";
		
		
		@Override
		protected void reduce(Text key, Iterable<Text> values,
				Context context)
				throws IOException, InterruptedException {
			
			//<Motif, Dataset_Targets>
			Map<String,SortedMap<String,Integer>> motifMap = new HashMap<String,SortedMap<String,Integer>>();
			
			
			
			for (Text t : values){
				
				String s = t.toString();
				String [] splits = s.split(paste);
				
				String motif = splits[0];
				String fileID= splits[1];
				int F = Integer.parseInt(splits[2]);
					
				
				SortedMap<String,Integer> matches = motifMap.get(motif);
				if (matches==null){
					matches = new TreeMap<String,Integer>();
					motifMap.put(motif,matches);
					matches.put(fileID,F);
				} else {
					Integer currentF = matches.get(fileID);
					
					if (currentF==null){
						matches.put(fileID,F);
					} else {
						matches.put(fileID,(int)Math.max(currentF,F));
					}
				}
			}
		
			for (Map.Entry<String,SortedMap<String,Integer>> e : motifMap.entrySet()){
				if (e.getValue().size()!=2){ //only motifblspairs in intersection
					continue;
				}
				
				SortedMap<String,Integer> familyOccValues = e.getValue();
				
				int first = familyOccValues.get(familyOccValues.firstKey());
				int second = familyOccValues.get(familyOccValues.lastKey());
				
				outputKey.set(e.getKey()+paste+key.toString().split(paste)[1]); //motif_BLS
				outputValue.set(second-first); //F2-F1
				
				context.write(outputKey,outputValue);
				
			}
			
			
			
			
			
			motifMap.clear();
			motifMap=null;
		}
	}
	
	
	
	
	@Override
	public int run(String[] args) throws Exception {
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "TargetComparator");
		
		job.setJarByClass(DBTargetComparator.class);
		
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
		job.setMapperClass(TargetMapper.class);
		job.setReducerClass(TargetReducer.class);

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
		System.out.println("bin/hadoop jar MaxBLSComparator be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DBMaxBLSComparator(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
	
	
}
