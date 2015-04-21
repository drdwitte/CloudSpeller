package be.iminds.cloudspeller.driver;

import be.iminds.cloudspeller.input.GeneFamilyInputFormat;

import java.io.IOException;
import java.util.Set;


import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFreqVec;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContentFactory;
import be.iminds.cloudspeller.motifpermutationgroups.MotifPermutationGroup;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.BytesWritable;

import org.apache.hadoop.io.Text;

import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.ToolRunner;

import be.iminds.cloudspeller.phylogenetics.BLS;

public class MotifDistributionAnalyzer extends frameWorkTest {
		
	public static class BLSDistribReducer extends Reducer<BytesWritable,BytesWritable, Text, Text> {

		protected final static int numberOfTries = 10000;
		protected static Text outputKey = new Text();
		protected static  Text outputValue = new Text();
 		private static CloudSpeller cloudSpeller;
 		//private static double [] alpha = { 0.001, 0.005, 0.01, 0.025, 0.05, 0.1 };

		@Override
		protected void reduce(BytesWritable inputKey, Iterable<BytesWritable> inputValues,
				Context context) throws IOException, InterruptedException {
			
			MotifContentFactory contentFact = cloudSpeller
					.getMotifContentFactory();
			MotifContent key = contentFact.createMotifContentFromBytes(inputKey);

			
			MotifPermutationGroup permGroup = new MotifPermutationGroup(key);
		
			// join frequency vectors belonging to the same motif and add to
			// motif map
			for (BytesWritable value : inputValues) {
				MotifFreqVec mFreq = MotifFreqVec.createMotifFreqVecFromBytes(value);
				permGroup.addMotifFreq(mFreq);
			}
			
			outputKey.set(key.toString());
			StringBuilder sb = new StringBuilder();
			
			//headerline
			//sb.append("\n (BLS \t mean \t sigma"); // \t alph=");
			/*for (double a : alpha){
				sb.append(a+"\t");
			}*/
			sb.append("\n");
			
			Set<Motif> background = key.createPermutationGroup(numberOfTries);
			
			double [] mean = permGroup.calculateMean(background);
			double [] mean2 = permGroup.calculateMeanSquare(background);
			double [] sigma = permGroup.calculateVariation(mean,mean2);
						
			for (int i=0; i<mean.length; i++){
				sb.append(BLS.getBLSThresholds()[i]+"\t");
				sb.append(mean[i]+"\t");
				sb.append(sigma[i]+"\t");
				/*int [] pvalOcc = permGroup.calculatePValueFamOcc(alpha,i,background);
				for (int j=0; j<pvalOcc.length; j++){
					sb.append(pvalOcc[j]+"\t");
				}*/
				sb.append("\n");
			}
			
			 
			
			
			
			//normality test
			/*sb.append("NormalityTest: BLS \t P_1sig=0.68 \t p_2sig=0.9545 p_3sig=0.9973 \n");
			
			
			double [] oneSig = permGroup.getPercentageInSigmaInterv(1,background,mean,sigma);
			double [] twoSig = permGroup.getPercentageInSigmaInterv(2,background,mean,sigma);
			double [] threeSig = permGroup.getPercentageInSigmaInterv(3,background,mean,sigma);
		
			for (int i=0; i<FreqVec.getNumberOfIntervals(); i++){	
				
				sb.append(BLS.getBLSThresholds()[i]+"\t"+oneSig[i]+"\t"+twoSig[i]+"\t"+threeSig[i]+"\n");
			}*/
			
			//nog een idee zijn percentielen maar dat moet in be.iminds.cloudspeller.postprocessing gebeuren!
			//we hebben sowieso al de top zoveel percent per permutatiegroep kunnen we 
			//ook een simulatie lopen met de top x % 
			
			
			outputValue.set(sb.toString());
			
			context.write(outputKey, outputValue);
			
		


			
		}


		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			cloudSpeller = new CloudSpeller(context.getConfiguration());
			cloudSpeller.initializeReducerFactories();
		}
 		
 		
 	}
	
	@Override
  	public int run(String[] args) throws Exception {
        
         //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = Job.getInstance(getConf(), "MotifDistributionAnalysis");
		
		job.setJarByClass(MotifDistributionAnalyzer.class);
		
   	  	//SETTING MAPREDUCE CLASSES
  	  	
		// IO
		job.setInputFormatClass(GeneFamilyInputFormat.class);
		//job.setInputFormatClass(MultiGFGFInputFormat.class); //TODO test op singlenodecluster werkt niet op amazon??
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(BytesWritable.class);
		job.setMapOutputValueClass(BytesWritable.class);

		// (K3,V3)
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(Text.class);

		// set hadoop methods
		job.setMapperClass(BLSMapper.class);
		//job.setCombinerClass(BLSCombiner.class);
		job.setReducerClass(BLSDistribReducer.class);

		processCommandLineArgs(args, job);
		
		printEnvironmentSettings(job.getConfiguration());
		
		
		
		//FIXME doesn't work on amazon
		
		// for multiple runs remove be.iminds.cloudspeller.output folder of previous run
		FileSystem fs = FileSystem.get(job.getConfiguration());
		// true stands for recursively deleting the folder you gave
		fs.delete(new Path(output), true);

		// setting IO

		try {
			FileInputFormat.setInputPaths(job, new Path(input));
			FileOutputFormat.setOutputPath(job, new Path(output));
		} catch (Exception e) {
			System.err.println(e.getMessage());
		}

		return job.waitForCompletion(true) ? 0 : 1;
  }
	
	





	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new MotifDistributionAnalyzer(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
	
	

}
