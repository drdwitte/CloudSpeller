package be.iminds.cloudspeller.postprocessing_Single;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import be.iminds.cloudspeller.motifmodels.FreqVec;


import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.SimultaneousOccurrenceConfidenceAndBLSFiltering;

import be.iminds.cloudspeller.phylogenetics.BLS;

public class DBMotifFilter extends Configured implements Tool {

	private static String input;
	private static String output;

	public static class MotifMapper extends Mapper<LongWritable,Text,NullWritable,Text> {

		private static final int [] t = {15,50,60,70,90,95};
				
		private static ConfidenceGraphRestrictions restrictions;
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			
			int localConfCutoff = context.getConfiguration().getInt("C",50);
			int localFamCutoff = context.getConfiguration().getInt("F",3);
			int localBLSCutoff = context.getConfiguration().getInt("BLS",0);
			//restrictions = new SimultaneousOccurrenceAndConfidenceFiltering(localConfCutoff,localFamCutoff);
			restrictions = new SimultaneousOccurrenceConfidenceAndBLSFiltering
				(localConfCutoff, localFamCutoff, localBLSCutoff );
		}

		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			String strValue = value.toString();
			if (strValue.length()==0)
				return;
			
			
			if (!restrictions.checkRestrictions(strValue)){
				return;
			}
			
			
				
			context.write(NullWritable.get(),value);

		}
		
	}
	
	@Override
	public int run(String[] args) throws Exception {
		
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "MotifFilter");
		
		job.setJarByClass(DBMotifFilter.class);
		
   	  	//SETTING MAPREDUCE CLASSES
  	  	
		// IO
		job.setInputFormatClass(TextInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(NullWritable.class);
		job.setMapOutputValueClass(Text.class);

		// (K3,V3)
		/*job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(IntWritable.class);*/

		// set hadoop methods
		job.setMapperClass(MotifMapper.class);
		job.setNumReduceTasks(0); //map only
		
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
		
		
		List<String> other_args = new ArrayList<String>();
		
		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-C".equals(args[i])) {
					confidenceCutoff = Integer.parseInt(args[++i]);
				} else if ("-F".equals(args[i])){
					familyCutoff = Integer.parseInt(args[++i]);
				
				} else if ("-BLS".equals(args[i])){
					minBLS = Integer.parseInt(args[++i]);	
					
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
		
		if (other_args.size()!=2){
			System.out.println("ERROR: Wrong number of parameters: " + args.length + " instead of 8.");
			printUsage();
		}
		
		input = other_args.get(0);
		output = other_args.get(1);
	}

	private static int printUsage() {
		System.out.println("bin/hadoop jar MotifFilter -C 50 -F 3 -BLS 0 be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DBMotifFilter(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
}
