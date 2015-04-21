package be.iminds.cloudspeller.driver;

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
import java.util.Set;

import be.iminds.cloudspeller.indexing.FullMotifMatch;
import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.NoDecorationFactory;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.indexing.PatternMatcher;
import be.iminds.cloudspeller.indexing.SequenceIDSet;
import be.iminds.cloudspeller.indexing.SuffixInterpreter;

import be.iminds.cloudspeller.indexing.PatternBLSPair;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.GeneFamilyInputFormat;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;
import be.iminds.cloudspeller.phylogenetics.BLSCalculatorFactory;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculatorFactory;
//FIXME nog eens kijken naar URIs want dit zal niet werken op amazon! (zie ook KN1Predictions of FrameworkTest)
public class DistributedABPatternMatcher extends Configured implements Tool {

	protected String input;
	protected String output;
	
	public static class GFABPatternMatcher extends Mapper<LongWritable, GeneFamily, Text,Text> {

		private static Set<PatternBLSPair> patterns;
		private static final int maxPatternLength = 15;

		private static Text outputKey = new Text();
		private static Text outputValue = new Text();

		private static ConservationScoreCalculatorFactory calculatorFac;
		
		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			patterns.clear();
			patterns=null;
		}

		@Override
		protected void map(LongWritable key, GeneFamily geneFamily, Context context)
		throws IOException, InterruptedException {
			
			
			
			if (geneFamily.getGenes().size()>15){
				System.out.println(geneFamily.getFamilyName()+" contains "+ geneFamily.getGenes().size()+ "genes -> ignored");
				return;
			} else {
				System.out.println("PM in: "+geneFamily.getFamilyName());
			}
			
			
			GSTFactory factory = new GSTFactory(maxPatternLength,true,new NoDecorationFactory());				
			
			
			PatternMatcher patternMatcher;
			patternMatcher = factory.createIndexStructure(geneFamily.getSequences());
			
			
			 
			BLSCalculator calculator = (BLSCalculator) calculatorFac.createCalculator(geneFamily);
			
			SuffixInterpreter interpreter = new SuffixInterpreter(geneFamily);
						
			patternLoop:for (PatternBLSPair p : patterns){
				
				List<Suffix> suffixes = patternMatcher.matchExactPattern(new IUPACMotif(p.getPattern()));
				BLS cutoffScore = new BLS(p.getBls());
								
				if (suffixes==null){
					continue patternLoop;
				}
				
				Map<Integer,List<Suffix>> alignedSuff = getAlignedSuffixesBySequencePosition(suffixes);
				
				NodeDecoration nodeInfo = new SequenceIDSet(geneFamily.getNumberOfGenes());
				
				Set<FullMotifMatch> uniqueMatches = new HashSet<FullMotifMatch>();
				
				for (Map.Entry<Integer,List<Suffix>> e : alignedSuff.entrySet()){
					
					nodeInfo.processSuffixes(e.getValue());
					BLS score = (BLS) calculator.calculateScore(nodeInfo);
					
					if (score == null || score.compareTo(cutoffScore)<0 ){ //score > BLSMIN /score >=patternBLS
						continue;
					}
										
					
					outputKey.set(p.toString());
					
					for (Suffix s : suffixes){
						FullMotifMatch match = interpreter.translateSuffixWithPos(p.getPattern(),s,p.getBls());
						uniqueMatches.add(match);
					}
					
				}
				
				for (FullMotifMatch m : uniqueMatches){
					m.writeToText(outputValue);
					context.write(outputKey,outputValue);
				}
			}
			
			

			
			
			
		}

		public static Map<Integer, List<Suffix>> getAlignedSuffixesBySequencePosition(List<Suffix> suffixes) {
			Map<Integer, List<Suffix>> suffixesAtPos = new HashMap<Integer, List<Suffix>>(); 
			for (Suffix s : suffixes){
				addToSuffixMap(s.getSequencePosition(),s, suffixesAtPos);
			}
			return suffixesAtPos;
			
		}
		
		public static void addToSuffixMap(int seqPos, Suffix s, Map<Integer,List<Suffix>> map){
			List<Suffix> suffixesAtPos = map.get(seqPos);
			if (suffixesAtPos == null){
				suffixesAtPos = new ArrayList<Suffix>();
				map.put(seqPos,suffixesAtPos);
			}
			suffixesAtPos.add(s);
		}

		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {

			patterns = new HashSet<PatternBLSPair>();
			
			Path[] uris = DistributedCache.getLocalCacheFiles(context
				    .getConfiguration());
			
			String patternsFile = uris[0].toString();
	    
			
		    BufferedReader in = new BufferedReader(new FileReader(patternsFile));
		    String patternBLS;
		    while ((patternBLS = in.readLine()) != null) {
		    	//System.out.println(patternBLS);
		    	patterns.add(new PatternBLSPair(patternBLS));
		    }
		    
		    in.close();
		    
		    calculatorFac = new BLSCalculatorFactory(new BLS(0));
			
		}
	}
	

	
	
	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new DistributedABPatternMatcher(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}
	
	
	@Override
	public int run(String[] args) throws Exception {
        
        //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = Job.getInstance(getConf(), "MotifDiscovery");
		
		job.setJarByClass(DistributedABPatternMatcher.class);
		
  	  	//SETTING MAPREDUCE CLASSES
 	  	
		// IO
		job.setInputFormatClass(GeneFamilyInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);

		// (K1,V1) map-> K2,V2) reduce-> (K3,V3)

		// (K2,V2)
		job.setMapOutputKeyClass(Text.class);
		job.setMapOutputValueClass(Text.class);

		// (K3,V3)
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(Text.class);

		// set hadoop methods
		job.setMapperClass(GFABPatternMatcher.class);
		// conf.setCombinerClass(Reduce.class);
		job.setReducerClass(Reducer.class); //identity reducer

		processCommandLineArgs(args, job);
		
		FileSystem fs = FileSystem.get(job.getConfiguration());
		// true stands for recursively deleting the folder you gave
		fs.delete(new Path(output), true);

		
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
				} else if ("-patterns".equals(args[i])){
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
		System.out.println("bin/hadoop jar TestMotifHadoop.jar -patterns patterns.txt be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
	}

}

