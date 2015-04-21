package be.iminds.cloudspeller.KN1Analysis;

import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Scanner;


import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;


import be.iminds.cloudspeller.output.BLSConfidenceGraph;

import be.iminds.cloudspeller.toolbox.GeneralToolbox;

public class KN1PredictionsFromDeNovo extends Configured implements Tool{
	
	private String input;
	private String output;
	
	public static class OutputScanner extends Mapper<LongWritable,Text, NullWritable,Text> {

		private static String regexp;
		private static String regexpComp;
		//private static MotifBLSRestrictions restrictions;
		
		
		

		@Override
		protected void map(LongWritable key, Text inputValue, Context context)
				throws IOException, InterruptedException {
			
			String line = inputValue.toString();
			if (line.length()==0){
				return;
			}
			
			//BLSConfidenceGraph graph = new BLSConfidenceGraph(line,motifFactory);
			
			Scanner scanner = GeneralToolbox.generateScanner(line);
			scanner.next(); 
			String motif = scanner.next();
			scanner.close();
			
			if (KN1Toolbox.motifMatchesWithRegExp(motif,regexp) || KN1Toolbox.motifMatchesWithRegExp(motif,regexpComp)){
			
				context.write(NullWritable.get(), inputValue);
			}			
		}

		
		@SuppressWarnings("unused")
		private boolean graphMatchesWithRegExp(BLSConfidenceGraph graph, String reg)  {
			return KN1Toolbox.motifMatchesWithRegExp(graph.getMotif().toString(),reg);
			
		}

		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
		
			regexp = "TGA..GA..GA.";
			regexpComp = ".TC..TC..TCA";
		}
		
		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			
		}
	
	}

	@Override
	public int run(String[] args) throws Exception {
        
		Job job = new Job(getConf(), "MotifExtrationFromOutput");
		
		job.setJarByClass(KN1PredictionsFromDeNovo.class);
		
  	  	//SETTING MAPREDUCE CLASSES
 	  	
		// IO
		//no be.iminds.cloudspeller.input format set -> TextInputFormat
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2) (motif, minBLS  #families ) 
		job.setMapOutputKeyClass(NullWritable.class);
		job.setMapOutputValueClass(Text.class);

		// set hadoop methods
		job.setMapperClass(OutputScanner.class);		
		//job.setReducerClass(Reducer.class); (sorting is not necessary -> omit)
		job.setNumReduceTasks(0);


		processCommandLineArgs(args, job);
		
		try {
			FileInputFormat.setInputPaths(job, new Path(input));
			FileOutputFormat.setOutputPath(job, new Path(output));
		} catch (Exception e) {
			System.err.println(e.getMessage());
		}
		
		//TODO remove dit werkt niet op amazon!
		// for multiple runs remove be.iminds.cloudspeller.output folder of previous run
		//FileSystem fs = FileSystem.get(getConf());
		// true stands for recursively deleting the folder you gave
		//fs.delete(new Path(be.iminds.cloudspeller.output), true);

		return job.waitForCompletion(true) ? 0 : 1;
	}

	private void processCommandLineArgs(String[] args, Job job) {

		// process command line arguments
		List<String> other_args = new ArrayList<String>();
		

		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-r".equals(args[i])) {
					job.setNumReduceTasks(Integer.parseInt(args[++i]));
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
	}



	static int printUsage() {
		System.out.println("bin/hadoop jar OutputSifter.jar be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
	}
	
	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new KN1PredictionsFromDeNovo(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}

	
	




}
