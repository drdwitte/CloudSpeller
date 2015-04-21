package be.iminds.cloudspeller.KN1Analysis;

import be.iminds.cloudspeller.input.Gene;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import be.iminds.cloudspeller.output.BLSConfidenceGraph;
import be.iminds.cloudspeller.output.MotifBLSRestrictions;
import be.iminds.cloudspeller.phylogenetics.BLS;

import be.iminds.cloudspeller.toolbox.LineIterator;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import be.iminds.cloudspeller.alphabets.CharacterIterator;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;


/**
 * Rework be.iminds.cloudspeller.output from KN1PredictionsDeNovo -> add pwm scores + sort + generate file for distributed
 * pattern matcher
 * @author ddewitte
 *
 */

public class DegMatchesProcessor {

	private Map<String,Double> bsScore = new HashMap<String,Double>();
	private MotifFactory factory = new IUPACFactory(IUPACType.TWOFOLDSANDN);
	private MotifBLSRestrictions restrictions;
	
	static {
		int [] t = {15,50,60,70,90,95};
		BLS.initializeBLSConstants(t);
		FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
	}
	
	public DegMatchesProcessor(MotifBLSRestrictions restrictions){
		this.restrictions=restrictions;
	}
	
	
	public static List<StringBuilder> createAllDegMatchesWithRegExp(String s, String alph) {
		
		List<StringBuilder> oldA = new ArrayList<StringBuilder>();
		List<StringBuilder> newA = new ArrayList<StringBuilder>();
		
		oldA.add(new StringBuilder(""));
		
		for (int i=0; i<s.length(); i++){
			char c = s.charAt(i);
			
			if (c == '.'){
				for (StringBuilder b : oldA){
					
					for (int j=0; j<alph.length(); j++){
						
						StringBuilder copy = new StringBuilder(b);
						copy.append(alph.charAt(j));
						newA.add(copy);
					}
				}
				oldA.clear();
				oldA.addAll(newA);
				
				newA.clear();
				//oldA = newA;
				
			} else {
				for (StringBuilder b : oldA){
					b.append(c);
				}
			}
		}
		
		return oldA;
	}
	
	public void initializeBSMap(String regex, PWM pwm){
		List<StringBuilder> bindingSites = createAllDegMatchesWithRegExp(regex,"ACGT");
		
		
		
		for (StringBuilder b : bindingSites){
			bsScore.put(b.toString(),pwm.calculateMatrixMatch(b.toString()));
		}
		
	}
	
	/**
	 * Rework the be.iminds.cloudspeller.output of KN1Predictions de Novo (sort + add PWm scores
	 * Line: TGAAAGAACGAC	15	3	100.0
	 * @param args
	 */
	public void reworkKN1Predictions(String inputfilename, String outputScoreWorst, PWM pwm) 
			throws IOException {
		
		String regex = KN1Toolbox.regex;
		LineIterator it = new LineIterator(inputfilename);
		
		
		SortedSet<ScorePair> scoreTableWorstBS = new TreeSet<ScorePair>();
		//SortedSet<ScorePair> scoreTableAverageBS = new TreeSet<ScorePair>();
		initializeBSMap(regex,pwm);
		
		
		
		while (it.hasNext()){
			
			String line = it.next();
			if (line.length()==0){
				continue;
			}
			
			BLSConfidenceGraph graph = new BLSConfidenceGraph(line, factory);
						
			boolean [] restr = restrictions.checkRestrictions(graph);
			
			int index=-1;
			for (int i=0; i<restr.length; i++){
				if (restr[i]){
					index=i;
					break;
				}
			}
			
			if (index<0){
				continue;
			}
			
			
			String motifStr = takeComplementIfNeeded(graph.getMotif().toString(), regex);
			int minBLS = BLS.getBLSThresholds()[index];
			int f = graph.getFreq(index);
			double conf = graph.getProbValue(index);
			
			String key = motifStr+"\t"+minBLS+"\t"+f+"\t"+conf;
			
			Double valueW = calculatescoreWorstBS(motifStr); 
			scoreTableWorstBS.add(new ScorePair(key,valueW));
			
		}
		
		BufferedWriter writer;
		if (outputScoreWorst!=null){
			writer = new BufferedWriter(new FileWriter(outputScoreWorst));
			
			
			
			for (ScorePair p : scoreTableWorstBS){
				
				//DEBUG
				if (p.getValue()<0.0)
					break;
				
				writer.write(p.getKey()+"\t"+p.getValue()+"\n");
								
			}
			writer.close();
		}
		
		/*if (outputScoreAvg!=null){
			writer = new BufferedWriter(new FileWriter(outputScoreAvg));
			for (ScorePair p : scoreTableAverageBS){
				writer.write(p.getKey()+"\t"+p.getValue()+"\n");
			}
			writer.close();
		}*/
	
	}
	
	/**
	 * 
	 * @param m een string matching with regex or its complement
	 * @param regex
	 * @return the complement of m if required to match with regex
	 */
	private String takeComplementIfNeeded(String m, String regex) {
		
		if (m.matches(regex)){
			return m;
		} else {
			return (new IUPACMotif(m,0)).getComplement().toString();
		}
		
	}

	public Double calculatescoreWorstBS(String motif) {
		
		List<StringBuilder> bsMatches = getMatchesWithMotif(motif);
		double bestScore = Double.MAX_VALUE;
		
		for (StringBuilder b : bsMatches){
			double currentScore = bsScore.get(b.toString());
			if (currentScore < bestScore){
				bestScore = currentScore;
			}
		}
		
		return new Double(bestScore);
		
	}
	
/*private Double calculatescoreAverageBS(String motif) {
		
		List<StringBuilder> bsMatches = getMatchesWithMotif(motif);
		double sum = 0.0; 
		for (StringBuilder b : bsMatches){
			sum+=bsScore.get(b.toString());
		}
		return new Double(sum/bsMatches.size());
		
	}*/
	
	

	private List<StringBuilder> getMatchesWithMotif(String motif) {
		IUPACAlphabet alph = new IUPACAlphabet(IUPACType.FULL);
				
		List<StringBuilder> oldA = new ArrayList<StringBuilder>();
		List<StringBuilder> newA = new ArrayList<StringBuilder>();
		
		oldA.add(new StringBuilder(""));
		
		for (int i=0; i<motif.length(); i++){
			char c = motif.charAt(i);
			
			if (alph.getNumberOfMatchingCharacters(c)>1){
				for (StringBuilder b : oldA){
					
					CharacterIterator it = alph.getMatchingCharactersIterator(c);
					
					while (it.hasNext()){
					
						StringBuilder copy = new StringBuilder(b);
						copy.append(it.next());
						newA.add(copy);
					}
				}
				oldA.clear();
				oldA.addAll(newA);
				
				newA.clear();
				//oldA = newA;
				
			} else {
				for (StringBuilder b : oldA){
					b.append(c);
				}
			}
		}
		
		return oldA;
		
		
		
	}

	public static void generateInputFileForDBPatternMatcher(String inputfilename, String outputFilename)
			throws IOException {
		
		LineIterator it = new LineIterator(inputfilename);
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilename));
		
		while (it.hasNext()){
			String line = it.next();
			Scanner scanner = new Scanner(line);
			writer.write(scanner.next()+"\t"+scanner.nextInt()+"\n");
			scanner.close();
		}
		
		writer.close();
		
	}
	
	public static Set<DeNovoRecord> selectDeNovoRecords(double pwmThresholdScore, String deNovoFile) throws IOException {
		Set<DeNovoRecord> records = new HashSet<DeNovoRecord>();
		LineIterator iterator = new LineIterator(deNovoFile);
		
		while (iterator.hasNext()){
			DeNovoRecord record = new DeNovoRecord(iterator.next());
			if (record.getPwmScore()>=pwmThresholdScore){
				records.add(record);
			}
		}
		return records;
	}
	
	/**
	 * Mais family
	 * @param records
	 * @param PMOutputFile
	 * @return
	 * @throws IOException
	 */
	
	public static Map<Gene,Set<RecS>> innerJoinDeNovoAndPMRecords(Set<DeNovoRecord> records, String PMOutputFile) throws IOException{
		
		Map<Gene,Set<RecS>> innerJoin = new HashMap<Gene,Set<RecS>>();
		
		Map<String,Double> filteredMotifs = new HashMap<String,Double>();
		
		//get all motif - score pairs 
		for (DeNovoRecord record : records){
			String motif = record.getMotif();
			filteredMotifs.put(motif,record.getPwmScore());
		}
				
		LineIterator iterator = new LineIterator(PMOutputFile);
		int lineCnt=0;
		while (iterator.hasNext()){
			lineCnt++;
			if (lineCnt%100000==0){
				System.out.print("*");
			}
			PMRecord record = new PMRecord(iterator.next());
			
			String geneStr = record.getGene();
			
			
			if (!geneStr.startsWith("ZM")){
				continue;
			}
			
			Gene gene = new Gene(geneStr,"ZM");
			String motif = record.getMotif();
			if (filteredMotifs.containsKey(motif)){
				
				Double score = filteredMotifs.get(motif);
				
				RecS smallRecord = new RecS(motif,score);
			
				Set<RecS> newRecords = innerJoin.get(gene);
				
				if (newRecords == null){
					newRecords = new HashSet<RecS>();
					innerJoin.put(gene,newRecords);
				}
				
				newRecords.add(smallRecord);
			}
			
			//DEBUG test
			//DEBUG test
			//DEBUG test
			//if (record.getMotif().equals("TGACTGACTGAC")){
			//	System.out.println(filteredMotifs.containsKey("Contains? "+ motif));
			//	System.out.println("Score: "+filteredMotifs.get(motif));
			//	
			//}
			//DEBUG test
			//DEBUG test
			//DEBUG test
			
		
		}
	
		return innerJoin;
	}

	
	
	public void generateInputFilesPM(String input, String outputRewWorstScore, String pmFile, PWM pwm) throws IOException{
		reworkKN1Predictions(input,outputRewWorstScore,pwm);
		generateInputFileForDBPatternMatcher(outputRewWorstScore, pmFile);
	}
	
	public static void reworkPMOutput(String outputRewWorstScore, String outputPM, String exp) throws IOException{
		int t = KN1Toolbox.getPWMScoreThreshold();
		String filename = KN1Toolbox.generateOutputFilenameDeNovoWorstBS(exp,t);
				
		Set<DeNovoRecord> kn1VariantsAbovePWMThreshold = selectDeNovoRecords(t,outputRewWorstScore);
		System.out.println("Number of De Novo records above PWM threshold " + t + " : "+kn1VariantsAbovePWMThreshold.size());
		
		Map<Gene,Set<RecS>> pmMatchesAboveThreshold = innerJoinDeNovoAndPMRecords(kn1VariantsAbovePWMThreshold,outputPM);
		System.out.println("Number of target maize genes matching with cumul above threshold : "+pmMatchesAboveThreshold.size());
		
		KN1Toolbox.printTableS(pmMatchesAboveThreshold,filename);
		
		
	}
	
	public static void main(String [] args) throws IOException {
		//KN1Toolbox.setKN1Directory("/home/ddewitte/Bureaublad/CloudSpellerExperimentsFinal/ABAnalysis/");
		
		
		String dir = KN1Toolbox.getKN1Directory();
		
		
		
		int [] F= {0};
		int [] C = {50,90};
		PWM pwm = KN1Toolbox.generateKN1PWM();
		
		//PWM pwm = KN1Toolbox.generateModifiedPWM();

		for (int f=0; f<F.length; f++){
			for (int c=0; c<C.length; c++){
						
			MotifBLSRestrictions restrictions = new MotifBLSRestrictions(F[f],C[c]);
			DegMatchesProcessor processor = new DegMatchesProcessor(restrictions);	
			String exp = "C"+C[c]+"F"+F[f]; // +"MODPWM";
		
			//String be.iminds.cloudspeller.input = dir +"DeNovoMatches6Nov13.txt";
			String input = dir  + "allABPMVariantMatches.txt";
			
			//String outputRewWorstScore = dir + "deNovoMatchesFiltered"+exp+".txt";
			String outputRewWorstScore = dir + "ABMatchesFiltered"+exp+".txt";
			//String pmFile = dir + "PatternsFiltered"+exp+".txt";
			String pmFile = dir + "ABPatternsFiltered"+exp+".txt";

			processor.generateInputFilesPM(input, outputRewWorstScore, pmFile,pwm);
			
			}
			
		}
		
		
		String outputPM = KN1Toolbox.getPMOutput(KN1Toolbox.getExperiment());
		
		String recordsFile = KN1Toolbox.getDeNovoRecords(KN1Toolbox.getExperiment());
		reworkPMOutput(recordsFile, outputPM,KN1Toolbox.getExperiment());

		
		
		
		
		
		
	}
	
	
	
	
		
		
		
		
		
	
	
	
	
	
	
}
