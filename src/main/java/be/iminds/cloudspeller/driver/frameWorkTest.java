package be.iminds.cloudspeller.driver;
//NOTE mapreduce libraries instead of mapred (deprecated)

import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.GeneFamilyInputFormat;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import java.util.List;

import be.iminds.cloudspeller.motifalgorithms.DiscoveryAlgorithm;
import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFreqVec;

import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContentFactory;
import be.iminds.cloudspeller.motifpermutationgroups.MotifPermutationGroup;
import be.iminds.cloudspeller.motifpermutationgroups.RedundanceChecker;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.BytesWritable;
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

import be.iminds.cloudspeller.output.ConfidenceGraph;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;

import be.iminds.cloudspeller.phylogenetics.ConservationScore;

public class frameWorkTest extends Configured implements Tool {

	protected URI uri;
	protected String input;
	protected String output;
	
	/**
	 * For each gene family the discovery algorithm is run. Every word satisfying the 
	 * constraints set in Configuration is emitted as:
	 * (Key,Value) = (PermutationGroupIdentifier,MotifFrequencyVector)
	 * @author Dieter De Witte
	 *
	 */
	public static class BLSMapper extends Mapper<LongWritable, GeneFamily, BytesWritable,BytesWritable> {

		private static MotifExtractor emitter;
		private static CloudSpeller cloudSpeller;
		private static int mapId=0;
		
		public void map(LongWritable key, GeneFamily geneFamily, Context context)
				throws IOException, InterruptedException {

			long startReading = System.nanoTime();
			if (geneFamily.getGenes().size()>15){
				System.out.println(geneFamily.getFamilyName()+" contains "+ geneFamily.getGenes().size()+ "genes -> ignored");
				return;
			} else {
				System.err.println("Mapping: "+geneFamily.getFamilyName()+ "\t"+(++mapId));
			}

			System.err.println("Reading (ms): "+(System.nanoTime()-startReading)/1000000);
			
			//TODO remove
			long start= System.nanoTime();
			HeapAnalyzer analyzer = new HeapAnalyzer();
			
			DiscoveryAlgorithm alg = cloudSpeller.getDiscoveryAlgorithm(geneFamily);

			alg.setMotifExtractor(emitter);
			alg.runDiscovery(cloudSpeller.getMotifFactory());
			
			analyzer.printAllNoMessage();
			long stop=System.nanoTime();
			System.err.println("Time (ms): "+(stop-start)/1000000);
		}
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			
			cloudSpeller = new CloudSpeller(context.getConfiguration());
			System.out.println(cloudSpeller);
			cloudSpeller.initializeMapperFactories();
			
			switch (cloudSpeller.getMotifAlgType()){
			
				case EXACTAB: 
					emitter = new ABEmitter(context);
					break;
				
				default:
					emitter = new InstantEmitter(context);
					break;
			}
			
			
		}
		
		@Override
		protected void cleanup(Context context) throws IOException,
				InterruptedException {
			
			emitter.close();
			System.err.println(emitter.getNumberOfMotifsExtracted()+ " motifs emitted");
		}
		
		
		
		
		class InstantEmitter implements MotifExtractor {
			
			private BytesWritable outputKey = new BytesWritable();
			private BytesWritable outputValue = new BytesWritable();
			private Context context;
			private int numberOfMotifsExtracted;
				
			private RedundanceChecker redundanceChecker = new RedundanceChecker(cloudSpeller.getMotifContentFactory(),cloudSpeller.getMotifFactory().getAlphabet());
			
			public InstantEmitter(Context context){
				this.context=context;
			}
			
			@Override
			public void add(Motif motif, ConservationScore score){
				
				numberOfMotifsExtracted++;
				MotifContentFactory contentFactory = cloudSpeller.getMotifContentFactory();
				// create permutation group ID
				MotifContent groupID = contentFactory.createContent(motif);

				if (redundanceChecker.isRedundant(groupID)){
					return;
				}
				
				// create motif frequency vector
				MotifFreqVec motifBLS = new MotifFreqVec(motif, score.createFrequencyVector());

				outputKey.set(contentFactory.createBytesRepresentation(groupID));
				outputValue.set(motifBLS.createBytesRepresentation());
				
				try {
					context.write(outputKey, outputValue);
				} catch (IOException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

			@Override
			public int getNumberOfMotifsExtracted() {
				return numberOfMotifsExtracted;
			}

			@Override
			public void reset() {
				
			}

			@Override
			public void close() {
				
			}


			
		}
		
		
		class PartialInMapperCombiner implements MotifExtractor {

			private BytesWritable outputKey = new BytesWritable();
			private BytesWritable outputValue = new BytesWritable();
			private Context context;
			private int numberOfMotifsExtracted=0;
			private int numberOfLocalCombinations=0;
			
			private int maxEquivalentLength = 7; 
			private Map<Motif,MotifFreqVec> shortMotifMap = new HashMap<Motif,MotifFreqVec>();
			
			private RedundanceChecker redundanceChecker = new RedundanceChecker(cloudSpeller.getMotifContentFactory(),cloudSpeller.getMotifFactory().getAlphabet());

			
			public PartialInMapperCombiner(Context context) {
				this.context = context;
			}

			@Override
			public void add(Motif motif, ConservationScore score) {
				int keq = motif.getGeneralizedLength();
				
				
				if (keq > maxEquivalentLength){ //instant emit
					 
					emitMotif(motif,new MotifFreqVec(motif,score.createFrequencyVector()));
					
				} else { //store
					MotifFreqVec mFreq =shortMotifMap.get(motif);
					FreqVec newFreqVec = score.createFrequencyVector();
					if (mFreq == null){
						shortMotifMap.put(motif,new MotifFreqVec(motif,newFreqVec));
					} else {
						numberOfLocalCombinations++;
						mFreq.getVec().add(newFreqVec);
					}
					
					if (shortMotifMap.size()%100000==0){
						System.out.println("Shortmotifmapsize: "+shortMotifMap.size());
					}
				}
			}

			@Override
			public void close() {
				for (Map.Entry<Motif,MotifFreqVec> e : shortMotifMap.entrySet()){

					emitMotif(e.getKey(),e.getValue());
				}
				
				System.out.println("number of motifs locally stored: "+shortMotifMap.size());
				shortMotifMap = null;
				
				System.out.println("number of motifs combined: "+numberOfLocalCombinations);
			}
			
			private void emitMotif(Motif m, MotifFreqVec mFreq){
					
				numberOfMotifsExtracted++;
				
				// create permutation group ID
				MotifContentFactory contentFactory = cloudSpeller.getMotifContentFactory();

				
				MotifContent groupID = contentFactory.createContent(m);
				
				if (redundanceChecker.isRedundant(groupID)){
					return;
				}
								
				outputKey.set(contentFactory.createBytesRepresentation(groupID));
				outputValue.set(mFreq.createBytesRepresentation());
								
				try {
					context.write(outputKey, outputValue);
				} catch (IOException ex) {
					ex.printStackTrace();
				} catch (InterruptedException ex) {
					ex.printStackTrace();
				}
			}

			@Override
			public int getNumberOfMotifsExtracted() {
				return numberOfMotifsExtracted;
			}

			@Override
			public void reset() {
				
			}


			
		}
		
		/**
		 * Only emits after a full sequence is scanned since a motif may have multiple
		 * occurrences in an alignment!
		 * @author ddewitte
		 *
		 */
		//FIXME uptdate ABMotifContainer!
		class ABEmitter implements ABMotifExtractor {

			private BytesWritable outputKey = new BytesWritable();
			private BytesWritable outputValue = new BytesWritable();
			
			private Context context;
			
			int numberOfMotifsExtracted=0;
			
			private Map<Motif,ConservationScore> motifMap = new HashMap<Motif,ConservationScore>();
			
			private RedundanceChecker redundanceChecker = new RedundanceChecker(cloudSpeller.getMotifContentFactory(),cloudSpeller.getMotifFactory().getAlphabet());



			public ABEmitter(Context context){
				this.context=context;
			}
			

			
			@Override
			/**
			 * NOTE: pay attention to reverse motif must also be added!!
			 */
			public void add(Motif motif, ConservationScore score) {
				
				ConservationScore oldScore = motifMap.get(motif);
				
				
				if (oldScore == null){
					motifMap.put(motif,score);
				}
				else if (oldScore.compareTo(score)<0){
					motifMap.put(motif,score);
					
				} else {}
				
			}

			@Override
			public void close() {
				//System.out.println("TempMapSize: "+motifMap.size());
				for (Map.Entry<Motif,ConservationScore> e : motifMap.entrySet()){
				
					emitMotif(e.getKey(),new MotifFreqVec(e.getKey(),e.getValue().createFrequencyVector()));			
					
					
				}
				
				motifMap.clear();
			}

			@Override
			public int getNumberOfMotifsExtracted() {
				return numberOfMotifsExtracted;
			}

			@Override
			/**
			 * Clear motif map en emit motifs (after maptask has completed
			 */
			public void reset() {
				
				this.close();
			}
			
			
			private void emitMotif(Motif m, MotifFreqVec mFreq){
				
				numberOfMotifsExtracted++;
				
				// create permutation group ID
				MotifContentFactory contentFactory = cloudSpeller.getMotifContentFactory();

				
				MotifContent groupID = contentFactory.createContent(m);
				
				if (redundanceChecker.isRedundant(groupID)){
					return;
				}
								
				outputKey.set(contentFactory.createBytesRepresentation(groupID));
				outputValue.set(mFreq.createBytesRepresentation());
								
				try {
					context.write(outputKey, outputValue);
				} catch (IOException ex) {
					ex.printStackTrace();
				} catch (InterruptedException ex) {
					ex.printStackTrace();
				}
			}
		
		}
		
		
		
	}
   
	
	public static class BLSCombiner extends Reducer<BytesWritable,BytesWritable,BytesWritable,BytesWritable> {
		
 		private static CloudSpeller cloudSpeller;

		private static BytesWritable outputValue = new BytesWritable();
		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
		
			cloudSpeller = new CloudSpeller(context.getConfiguration());
			cloudSpeller.initializeCombinerFactories();
		}
		
		public void reduce(BytesWritable inputKey, Iterable<BytesWritable> byteValues,
				Context context) throws IOException, InterruptedException {
		
			
			MotifContentFactory contentFact = cloudSpeller.getMotifContentFactory();
			MotifContent key = contentFact.createMotifContentFromBytes(inputKey);
	
			MotifPermutationGroup permGroup = new MotifPermutationGroup(key);

			// join frequency vectors belonging to the same motif and add to
			// motif map
			int numberIn=0;
			MotifFreqVec mFreq = new MotifFreqVec();
			for (BytesWritable value : byteValues) {
				mFreq = MotifFreqVec.createMotifFreqVecFromBytes(value);
				permGroup.addMotifFreq(mFreq);
				numberIn++;
			}
			
			//emit all frequency vectors
			int numberOut=0;
			for (Map.Entry<Motif,FreqVec> mf : permGroup.getAggregatedMotifMap().entrySet()){
							
				
				mFreq.setMotif(mf.getKey()); //recycle motifFrequency vectory
				mFreq.setVec(mf.getValue());
				
				outputValue.set(mFreq.createBytesRepresentation());
				
				context.write(inputKey, outputValue);
				numberOut++;
			}
			
			System.err.println("i/o="+numberIn + "/" + numberOut);
			
			
		}
		
		
	}
	
	/**
	 * Reducer processes all motifFrequency vectors belonging to a single 
	 * permutation group. For this group a background model is built and all 
	 * motifs satisfying the confidence constraints are written to the be.iminds.cloudspeller.output
	 * @author Dieter De Witte
	 *
	 */
 	public static class BLSReducer extends Reducer<BytesWritable,BytesWritable, Text, Text> {
		
 		protected static Text outputKey = new Text();
		protected static  Text outputValue = new Text();
 		private static CloudSpeller cloudSpeller;
		
 		
		@Override
		protected void setup(Context context) throws IOException,
				InterruptedException {
			
			cloudSpeller = new CloudSpeller(context.getConfiguration());
			cloudSpeller.initializeReducerFactories();
		}
  
		public void reduce(BytesWritable byteKey, Iterable<BytesWritable> byteValues,
				Context context) throws IOException, InterruptedException {

			long start = System.nanoTime();
			
			//TODO remove
			HeapAnalyzer analyzer = new HeapAnalyzer();
						
			MotifContentFactory contentFact = cloudSpeller
					.getMotifContentFactory();
			MotifContent key = contentFact.createMotifContentFromBytes(byteKey);

			
			MotifPermutationGroup permGroup = new MotifPermutationGroup(key);
		
			// join frequency vectors belonging to the same motif and add to
			// motif map
			for (BytesWritable value : byteValues) {
				MotifFreqVec mFreq = MotifFreqVec.createMotifFreqVecFromBytes(value);
				permGroup.addMotifFreq(mFreq);
			}
			
			//TODO remove
			System.err.println(key.toString()+" unique be.iminds.cloudspeller.input:" + "\t" +permGroup.groupSize());
			analyzer.printUsedMemory(); 

			// calc background model for permutation group
			permGroup.generateBackgroundModel();

			ConfidenceGraphRestrictions restrictions = cloudSpeller.getConfGraphRestrictions();
			OutputExtractor extractor = new ConfGraphExtractor(context, restrictions, key);
			extractor.setKeyValueTextContainers(outputKey, outputValue);
			
			permGroup.generateConfidenceGraphs(cloudSpeller.getConfGraphFactory(),extractor);

			//TODO remove
			System.err.println("Reduce time (ms): " + (System.nanoTime()-start)/1000000);
			
		}
		
		class ConfGraphExtractor implements OutputExtractor {
			
			private Text outputKey; 
			private Text outputValue;
			
			private int numMotifsExtracted=0;
			
			
			private ConfidenceGraphRestrictions restrictions;
			private Context context;
			private MotifContent key;
			
			public ConfGraphExtractor(Context context, ConfidenceGraphRestrictions restrictions, MotifContent key){
				this.restrictions = restrictions;
				this.context = context;
				this.key = key;
			}
			
			public void setKeyValueTextContainers(Text outputKey, Text outputValue){
				this.outputKey=outputKey;
				this.outputValue=outputValue;
					
			}
			
			public void extract(ConfidenceGraph graph) {
				
				if (++numMotifsExtracted%100000 == 0){
					System.err.println(numMotifsExtracted+" motifs extracted!");
				}
				try {
					if (restrictions.checkRestrictions(graph)) {
						outputKey.set(key.toString());
						graph.toText(outputValue);
						context.write(outputKey, outputValue);
						
					}
				} catch (Exception e){
					System.err.println("Exception thrown while writing graph for key: "+key.toString());
					e.printStackTrace();
				}
			}
		}
		
		class IgnoreOutput implements OutputExtractor {
		
			public IgnoreOutput(){
			}
			
			public void setKeyValueTextContainers(Text outputKey, Text outputValue){
			}
			
			public void extract(ConfidenceGraph graph) {
			}
		}
	}
  
	protected void processCommandLineArgs(String[] args, Job job) {

		// process command line arguments
		List<String> other_args = new ArrayList<String>();
		
		String settingsFile=null;
		
		for (int i = 0; i < args.length; ++i) {
			try {
				if ("-r".equals(args[i])) {
					job.setNumReduceTasks(Integer.parseInt(args[++i]));
				} else if ("-settings".equals(args[i])){
					settingsFile=args[++i];
				} else if ("-uri".equals(args[i])){
					uri=new URI(args[++i]);
				} else {
					other_args.add(args[i]);
				}
			} catch (NumberFormatException except) {
				System.out.println("ERROR: Integer expected instead of " + args[i]);
				printUsage();
			} catch (ArrayIndexOutOfBoundsException except) {
				System.out.println("ERROR: Required parameter missing from " + args[i - 1]);
				printUsage();
			} catch (URISyntaxException e) {
				System.out.println("URI exception");
				e.printStackTrace();
			}
		}
		
		if (settingsFile!=null){
			processSettingsFile(uri, settingsFile,job);
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
  
	//FIXME depends on deprecated function FSDataInputStream.readLine()!!
	@SuppressWarnings("deprecation")
	protected void processSettingsFile(URI uri, String file, Job job) {
		Configuration conf = job.getConfiguration();
        
		
		Path inFile = new Path(file);
		
		try {
			FileSystem fs;
			if (uri != null) {
				fs = FileSystem.get(uri, conf);
			} else {
				fs = FileSystem.get(conf);
			}
			if (!fs.exists(inFile))
				System.out.println("Unable to open settings file: "+file);
			
			FSDataInputStream in = fs.open(inFile);
			//FIXME  -> still needs to be tested, remove previous line
			//BufferedReader in = new BufferedReader(new InputStreamReader(fs.open(inFile)));
			String line;
			
			while ((line = in.readLine())!=null){
				
				String [] splits = line.split("=");
				if (line.startsWith("s")){
					conf.set(splits[0].substring(1),splits[1]);
					
				} else if (line.startsWith("i")){
					conf.setInt(splits[0].substring(1),Integer.parseInt(splits[1]));
				
				} else if (line.startsWith("d")){
					conf.setFloat(splits[0].substring(1),Float.parseFloat(splits[1]));
				}
	
			}
			
			in.close();
			
		} catch (IOException e) {
			System.out.println("Exception thrown while processing settingsfile");
			e.printStackTrace();
		}
		
		
		

	}

	public static void main(String[] args) throws Exception {

		int res = 0;
		try {
			res = ToolRunner
					.run(new Configuration(), new frameWorkTest(), args);
		} catch (Exception ex) {
			System.err.println("exception: " + ex.getMessage());
		}
		System.exit(res);
	}

  @Override
  	public int run(String[] args) throws Exception {
        
         //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
		Job job = Job.getInstance(getConf(), "MotifDiscovery");
		
		job.setJarByClass(frameWorkTest.class);
		
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
		job.setReducerClass(BLSReducer.class);

		processCommandLineArgs(args, job);
		
		printEnvironmentSettings(job.getConfiguration());
			
		
		//FIXME doesn't work on amazon
		
		// for multiple runs remove be.iminds.cloudspeller.output folder of previous run
		//FileSystem fs = FileSystem.get(job.getConfiguration());
		// true stands for recursively deleting the folder you gave
		//fs.delete(new Path(be.iminds.cloudspeller.output), true);

		// setting IO

		try {
			FileInputFormat.setInputPaths(job, new Path(input));
			FileOutputFormat.setOutputPath(job, new Path(output));
		} catch (Exception e) {
			System.err.println(e.getMessage());
		}
		
		System.err.println("Current Number of reduce tasks: "+job.getNumReduceTasks());

		//TODO remove
		/*if (job.getNumReduceTasks()==1){
			job.setNumReduceTasks(60);
			System.err.println(job.getNumReduceTasks());
		}*/
		
		return job.waitForCompletion(true) ? 0 : 1;
  }

	protected void printEnvironmentSettings(Configuration configuration) {
		
		//String notSetStr = "NOT_SET";
		//int notSetInt = 0;
		String MB = " MB";
		
		int memoryDataNode = 1000;
		int memoryTaskTracker = 1000;
		
		int numberOfMapSlots=configuration.getInt("mapred.tasktracker.map.tasks.maximum",0);
		int numberOfReduceSlots=configuration.getInt("mapred.tasktracker.reduce.tasks.maximum",0);
			
		System.err.println("************************************************************");
		
		String memMapSlot = configuration.get("mapred.map.child.java.opts","-Xmx200m");
		String memRedSlot = configuration.get("mapred.reduce.child.java.opts","-Xmx200m");
		
		int mbInMapSlot = parseVM(memMapSlot);
		int mbInRedSlot = parseVM(memRedSlot);
		
		
		int total = memoryDataNode + memoryTaskTracker;
		total+=numberOfMapSlots*mbInMapSlot;
		total+=numberOfReduceSlots*mbInRedSlot;
		
		System.err.println("************************************************************");
		System.err.println("Memory per datanode (Hadoop ref p 305): ");
		System.err.println("************************************************************");
		System.err.println("Datanode: "+memoryDataNode+MB);
		System.err.println("TaskTracker: "+memoryTaskTracker+MB);
		System.err.println("Tasktracker child map task (#slots x mem): "+numberOfMapSlots+" x "+memMapSlot+MB);
		System.err.println("Tasktracker child reduce task (#slots x mem): "+numberOfReduceSlots+" x "+memRedSlot+MB);
		System.err.println("Total: "+total+MB);
		

	
	}

	private int parseVM(String s) {
		if (s.length()==0){
			return 0;
		}
		
		int numMbPerGb=1024;
		int numMb = 0;
		String [] arguments = s.split(",");
		for (String arg : arguments){
			
			if (arg.startsWith("-Xmx")){
				numMb = Integer.parseInt(arg.substring(4,arg.length()-1));
				if (arg.charAt(arg.length()-1)=='m'){
					
				} else {//Gigabytes! 
					numMb*=numMbPerGb;
				}
			}
				
			
		}
		return numMb;
	}

	protected static int printUsage() {
		System.out.println("bin/hadoop jar TestMotifHadoop.jar -settings settingsfile be.iminds.cloudspeller.input/ be.iminds.cloudspeller.output/");
		return -1;
	}
	
	
	
	
  
}