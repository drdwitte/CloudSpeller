package be.iminds.cloudspeller.postprocessing_Single;


import java.io.IOException;
import java.util.Scanner;
import java.util.SortedSet;
import java.util.TreeSet;

import be.iminds.cloudspeller.motifmodels.FreqVec;

import org.apache.hadoop.conf.Configured;
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

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import be.iminds.cloudspeller.KN1Analysis.ScorePair;


public class DBMotifScoreSorter extends Configured implements Tool {

	private static String input;
	private static String output;
	
	
	static class MotifScoreSorterMapper extends Mapper<LongWritable,Text,DoubleWritable,Text> {
	
		//private static MotifFactory motifFactory;
		private static final int [] t = {15,50,60,80,90,95};
		
		private static DoubleWritable outputKey = new DoubleWritable();
		private static Text outputValue = new Text();
		
		private static  double maxCat=10.0;
		private static double step =0.2;
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			BLS.initializeBLSConstants(t);
			FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
			//motifFactory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
			
		}
	
	
		@Override
		protected void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			
			
			Scanner scanner = GeneralToolbox.generateScanner(value.toString());
			
			String motifBLS = scanner.next();
			double score = scanner.nextDouble();
			scanner.close();
			
			emitInCategory(motifBLS,score, context);
			
			
		}
	
		private static void emitInCategory(String motifBLS, double score, Context context) throws IOException, InterruptedException {
			
			outputValue.set(motifBLS+"\t"+score);
			
			
			double category=0.0;
			loop:for (double cat=1.0; cat<maxCat; cat+=step){
				if (score<cat){
					category=cat;
					break loop;
				}
			}
			if (category<1.0){
				category=maxCat;
			}
			
			outputKey.set(category);
			context.write(outputKey,outputValue);
			
		}
	}
	
	static class MotifScoreSorterReducer extends Reducer<DoubleWritable,Text,Text,DoubleWritable> {
	
		private static Text outputKey = new Text();
		private static DoubleWritable outputValue = new DoubleWritable();

		@Override
		protected void reduce(DoubleWritable key, Iterable<Text> values,
				Context context) throws IOException, InterruptedException
		{
			
			SortedSet<ScorePair> highscores = new TreeSet<ScorePair>();
			
			for (Text t : values){
			
				String v = t.toString();
				Scanner scanner = GeneralToolbox.generateScanner(v);
				
				String motifBLS = scanner.next();
				double score = scanner.nextDouble();
						
				highscores.add(new ScorePair(motifBLS,score));
				
			}
			
			for (ScorePair p : highscores){ //sort them per reducer outputfile
				
				outputKey.set(p.getKey());
				outputValue.set(p.getValue());
				
				
				context.write(outputKey,outputValue);
				
			}
			
			
		}
	}
		
	

	@Override
	public int run(String[] args) throws Exception {
		
		//NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = new Job(getConf(), "MotifScoreSorter");
		
		job.setJarByClass(DBMotifScoreSorter.class);
		
   	  	//SETTING MAPREDUCE CLASSES
  	  	
		// IO
		job.setInputFormatClass(TextInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(DoubleWritable.class);
		job.setMapOutputValueClass(Text.class);

		// (K3,V3)
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(DoubleWritable.class);

		// set hadoop methods
		job.setMapperClass(MotifScoreSorterMapper.class);
		job.setReducerClass(MotifScoreSorterReducer.class);

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
		
		if (args.length!=2){
			System.out.println("ERROR: Wrong number of parameters: " + args.length + " instead of 2.");
			printUsage();
		}
		input = args[0];
		output = args[1];
	}

	private static int printUsage() {
		System.out.println("bin/hadoop jar SecondarySorter.jar be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
		
	}



	
}
