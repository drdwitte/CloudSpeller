package be.iminds.cloudspeller.KN1Analysis;

import be.iminds.cloudspeller.indexing.GeneralizedSuffixTree;
import be.iminds.cloudspeller.indexing.NoDecorationFactory;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.Sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import be.iminds.cloudspeller.toolbox.LineIterator;


import be.iminds.cloudspeller.motifmodels.IUPACMotif;

public class KN1Toolbox {

	public static final boolean withrevcomp = true;
	public static final int maxPatternLength = 12;
	public static final int maxGenesInFamily = 15;
	public static Map<String,GeneFamily> families;
	public static String referenceMotif = "TGANNGANNGAN";
	public static String regex = "TGA..GA..GA.";
	private static String KN1Directory = "/home/ddewitte/Bureaublad/CloudSpellerExperimentsFinal/KN1_6NOV/";
	private static String datasetfilename="allFamilies.txt";
	private static String experiment = "C90F0";
//
	private static int pwmThreshold=0;
//
	public static boolean motifMatchesWithRegExp(String motif, String regex){
		
		if (motif.length() == regex.length()){
			
			for (int i=0; i<regex.length(); i++){
				if (regex.charAt(i)=='.'){
					continue;
				} else if (regex.charAt(i) == motif.charAt(i)){
					continue;
				} else {
					return false;
				}
			}
			return true;
		} 
		else 
			return false;
	}
	
	public static void setKN1Directory(String d){
		KN1Directory=d;
	}
	
	public static String getPMOutput(String exp){
		return KN1Directory + "outputPMs/outputPM" + exp + "/part-r-00000";
	}
	
	public static String getDeNovoRecords(String exp){
		return KN1Directory + "DeNovoMatchesFiltered/deNovoMatchesFiltered"+exp+".txt";
		
	}
	
	public static String getKN1Directory(){
		return KN1Directory;
	}
	
	public static String getDatasetPath(){
		
		return KN1Directory+datasetfilename;
	}
	
	public static void setPWMScoreThreshold(int t){
		pwmThreshold = t;
	}
	
	public static int getPWMScoreThreshold() {
		
		return pwmThreshold;
	}
	


	private static String generateFilename(String theme, int t){
		String filename;
		if (t<0){
			filename = KN1Directory + theme + "m" + "Pwm" +  (-t) + ".txt";
		} else {
			filename = KN1Directory + theme + "Pwm"+ t + ".txt";;
		}
		
		return filename;
	}
	
	public static String generateOutputFilenameBolduc(int t){
		return generateFilename("MaizeGeneExpPWMMatch",t);
//		return generateFilename("MaizeGeneExpPWMMatchMODPWM",t);

	}
	
	public static String generateOutputFilenameBolducWithOrthologs(int t){
		return generateFilename("MaizeGenesExpWithPWMMatchInOrthologs",t);
//		return generateFilename("MaizeGenesExpWithPWMMatchInOrthologsMODPWM",t);

	}
	
	public static String generateOutputFilenameDeNovoWorstBS(String exp, int t){
		return generateFilename("MaizeGenesDeNovoWorstBS"+exp,t);
	}
	
	public static String generateOutputFilenamePWMScan(int t){
		return generateFilename("MaizeGEnesPWMScan",t);
//		return generateFilename("MaizeGEnesPWMScanMODPWM",t);

	}
	
	public static Map<String,GeneFamily> getFamilies() {
		return families;
	}

	public static void setDatasetFromFile(String filenameDataset) throws IOException {
		if (families != null){
			return;
		}
		
		families = new HashMap<String,GeneFamily>();
		BufferedReader in = new BufferedReader(new FileReader(filenameDataset));

		GeneFamily gf = new GeneFamily(in);
		while (gf.isInitialized()){
			
			if (gf.getGenes().size()<=maxGenesInFamily){
				families.put(gf.getFamilyName(),gf);
			}
			gf = new GeneFamily(in);
		}
		in.close();
	}
	
	public static String getActualBS(Suffix s, Sequence seq) {
		if (s.getSequenceID()==0){
			return seq.toString().substring(s.getSequencePosition(),s.getSequencePosition()+KN1Toolbox.maxPatternLength);
		} else { //rev strand
			int index = seq.length()-1-s.getSequencePosition();
			String bs = seq.toString().substring(index+1-maxPatternLength, index+1);
			return getComplementOf(bs);
		}
	}

	public static String getComplementOf(String bs) {
		StringBuilder sb = new StringBuilder();
		
		for (int i=bs.length()-1; i>=0; i--){
			sb.append(getComplementChar(bs.charAt(i)));
		}
		
		return sb.toString();
	
	}

	public static char getComplementChar(char c) {
		if (c=='A')
			return 'T';
		else if (c=='C'){
			return 'G';
		} else if (c=='G'){
			return 'C';
		} else {
			return 'A';
		}
}

	public static List<Suffix> giveMatchesOfPatternInSequence(String pattern, Sequence seq){
		IUPACMotif motif = new IUPACMotif(pattern);
		ArrayList<Sequence> seqs = new ArrayList<Sequence>();
		seqs.add(seq);
		GeneralizedSuffixTree gst = new GeneralizedSuffixTree(seqs, 
				withrevcomp, maxPatternLength, new NoDecorationFactory());
		return gst.matchExactPattern(motif);
	}
	
	
	public static PWM generateModifiedPWM(){
	
		PWM pwm = new PWM(12);
		
		double [] backgr = {0.275, 0.225, 0.225, 0.275};
		pwm.setBackgroundModel(backgr);
		
		double [] probs0 = {0.000000, 0.000000, 0.000000, 1.000000};
		pwm.setProbabilities(0, probs0);
		
		double [] probs1 = {0.000000, 0.000000, 1.000000, 0.000000};
		pwm.setProbabilities(1, probs1);
		
		double [] probs2 = {1.000000, 0.000000, 0.000000, 0.000000};
		pwm.setProbabilities(2, probs2);
		
		double [] probs3 = {0.000000, 0.820000, 0.000000, 0.180000};
		pwm.setProbabilities(3, probs3);
		
		double [] probs4 = {0.100000, 0.080000, 0.320000, 0.500000};
		pwm.setProbabilities(4, probs4);
		
		double [] probs5 = {0.000000, 0.000000, 1.000000, 0.000000};
		pwm.setProbabilities(5, probs5);
		
		double [] probs6 = {1.000000, 0.000000, 0.000000, 0.000000};
		pwm.setProbabilities(6, probs6);
		
		double [] probs7 = {0.000000, 0.380000, 0.000000, 0.620000};
		pwm.setProbabilities(7, probs7);
		
		double [] probs8 = {0.100000, 0.080000, 0.320000, 0.500000};
		pwm.setProbabilities(8, probs8);
		
		double [] probs9 = {0.000000, 0.000000, 1.000000, 0.000000};
		pwm.setProbabilities(9, probs9);
		
		double [] probs10 = {1.000000, 0.000000, 0.000000, 0.000000};
		pwm.setProbabilities(10, probs10);
		
		double [] probs11 = {0.000000, 0.820000, 0.000000, 0.180000 };
		pwm.setProbabilities(11, probs11);
		
		return pwm;
	}
	
	public static PWM generateKN1PWM(){
	
		//return generateModifiedPWM();
		
		PWM pwm = new PWM(12);
		
		double [] backgr = {0.275, 0.225, 0.225, 0.275};
		pwm.setBackgroundModel(backgr);
		
		double [] probs0 = {0.000000, 0.000000, 0.000000, 1.000000};
		pwm.setProbabilities(0, probs0);
		
		double [] probs1 = {0.000000, 0.000000, 1.000000, 0.000000};
		pwm.setProbabilities(1, probs1);
		
		double [] probs2 = {1.000000, 0.000000, 0.000000, 0.000000};
		pwm.setProbabilities(2, probs2);
		
		double [] probs3 = {0.000000, 0.820000, 0.000000, 0.180000};
		pwm.setProbabilities(3, probs3);
		
		double [] probs4 = {0.100000, 0.080000, 0.320000, 0.500000};
		pwm.setProbabilities(4, probs4);
		
		double [] probs5 = {0.000000, 0.000000, 1.000000, 0.000000};
		pwm.setProbabilities(5, probs5);
		
		double [] probs6 = {1.000000, 0.000000, 0.000000, 0.000000};
		pwm.setProbabilities(6, probs6);
		
		double [] probs7 = {0.000000, 0.380000, 0.000000, 0.620000};
		pwm.setProbabilities(7, probs7);
		
		double [] probs8 = {0.250000, 0.000000, 0.250000, 0.500000};
		pwm.setProbabilities(8, probs8);
		
		double [] probs9 = {0.000000, 0.000000, 1.000000, 0.000000};
		pwm.setProbabilities(9, probs9);
		
		double [] probs10 = {1.000000, 0.000000, 0.000000, 0.000000};
		pwm.setProbabilities(10, probs10);
		
		double [] probs11 = {0.000000, 0.820000, 0.000000, 0.180000 };
		pwm.setProbabilities(11, probs11);
		
		return pwm;
	}
	
	/**
	 * Read table (small records) from file
	 * @throws IOException 
	 */
	public static Map<Gene, Set<RecS>> readTableSFromFile(String filename) throws IOException{
		Map<Gene, Set<RecS>> table = new HashMap<Gene, Set<RecS>>();
		LineIterator lineIt = new LineIterator(filename);
		int c=0;
		while (lineIt.hasNext()){
			c++;
			if (c%100000==0){
				System.out.print("*");
			}
			String line = lineIt.next();
			Scanner scanner = new Scanner (line);
			scanner.useLocale(Locale.US);
			
			Gene queryGene = new Gene(scanner.next(),"ZM");
			String bindingSite = scanner.next();
			double pwmScore = scanner.nextDouble();
			
			scanner.close();
			
			Set<RecS> records = table.get(queryGene);
			if (records == null){
				records = new HashSet<RecS>();
				table.put(queryGene,records);
			}
			records.add(new RecS(bindingSite, pwmScore));
			
		}
		return table;
		
	}
	
	/**
	 * Read table (full records) from file
	 * @throws IOException 
	 */
	public static Map<Gene, Set<RecF>> readTableFFromFile(String filename) throws IOException{
		Map<Gene, Set<RecF>> table = new HashMap<Gene, Set<RecF>>();
		LineIterator lineIt = new LineIterator(filename);
		while (lineIt.hasNext()){
			String line = lineIt.next();
			Scanner scanner = new Scanner (line);
			scanner.useLocale(Locale.US);
			
			Gene queryGene = new Gene(scanner.next(),"ZM");
			String geneFamily = scanner.next();
			String gene = scanner.next();
			
			String bindingSite = scanner.next();
			double pwmScore = scanner.nextDouble();
			
			scanner.close();
			
			Set<RecF> records = table.get(queryGene);
			if (records == null){
				records = new HashSet<RecF>();
				table.put(queryGene,records);
			}
			records.add(new RecF(geneFamily, gene, bindingSite, pwmScore));
			
		}
		return table;
	}
	
	/**
	 * Add pwm scores to table
	 * @param tableSmall, tableFull if a table is not availaible set to null
	 */
	public static void generateWeightForAllMatches(Map<Gene, Set<RecS>> tableSmall, Map<Gene, Set<RecF>> tableFull, PWM pwm) {
	
		Map<String,Double> bsScoreMap = new HashMap<String,Double>();
		
		if (tableSmall!=null){
			for (Map.Entry<Gene, Set<RecS>> e : tableSmall.entrySet()){
				for (RecS record : e.getValue()){
					String bs = record.getBindingSite();
					Double score = bsScoreMap.get(bs);
					if (score!=null){
						record.setPWMScore(score);
					} else {
						record.setPWMScore(pwm.calculateMatrixMatch(bs));
					}
				}
			}
		}
		
		if (tableFull!=null){
			for (Map.Entry<Gene, Set<RecF>> e : tableFull.entrySet()){
				for (RecF record : e.getValue()){
					String bs = record.getBindingSite();
					Double score = bsScoreMap.get(bs);
					if (score!=null){
						record.setPWMScore(score);
					} else {
						record.setPWMScore(pwm.calculateMatrixMatch(bs));
					}
				}
			}
		}
	

	}
	
	public static Map<Gene, Set<RecS>> selectRecSAbovePWMThreshold(Map<Gene, Set<RecS>> tableSmall, double pwmThreshold){
		
		Map<Gene, Set<RecS>> filteredTable = new HashMap<Gene, Set<RecS>>();
		
		for (Map.Entry<Gene,Set<RecS>> e : tableSmall.entrySet()){
			
			Set<RecS> records = new HashSet<RecS>();
			
			for (RecS record : e.getValue()){
				if (record.getPWMScore() >=pwmThreshold){
					records.add(record);
				}
			}
			
			if (records.size()>0){
				filteredTable.put(e.getKey(),records);
			}
		}
		
		return filteredTable;
		
	}
	
	/**
	 * Print table: Query gene -  Binding site - Score
	 * @param filename
	 */
	public static void printTableS(Map<Gene, Set<RecS>> tableS, String filename) throws IOException {
		
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		
		for (Map.Entry<Gene, Set<RecS>> e : tableS.entrySet()){
			String queryGene = e.getKey().getID();
			System.out.print("*");
			for (RecS record : e.getValue()){
				out.write(queryGene);
				out.write("\t");
				out.write(record.toString());
				out.write("\n");
			}
		}
		System.out.println("");

		out.close();
		
	}
	
	public static void printTableF(Map<Gene, Set<RecF>> tableF, String filename) throws IOException {
		
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		
		for (Map.Entry<Gene, Set<RecF>> e : tableF.entrySet()){
			String queryGene = e.getKey().getID();
			System.out.print("*");
		
			for (RecS record : e.getValue()){
				out.write(queryGene);
				out.write("\t");
				out.write(record.toString());
				out.write("\n");
			}
			
		}
		System.out.println("");
		out.close();
	}

	public static String getExperiment() {
		return experiment;
	}
	
	public static void setExperiment(String exp){
		experiment=exp;
	}



	
	
	
	
}
