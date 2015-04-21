package be.iminds.cloudspeller.IPA1Analysis;

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
import be.iminds.cloudspeller.KN1Analysis.RecF;
import be.iminds.cloudspeller.KN1Analysis.RecS;

public class IPA1Toolbox {

	private static final int maxGenesInFamily=15;
	public static final String reference = "TGGGCY";
	public static final String species = "OS";
	private static final int maxPatternLength=reference.length();
	private static final boolean withrevcomp=true;
	
	private static String IPA1Directory = "/home/ddewitte/Desktop/CloudSpellerExperimentsFinal/be.iminds.cloudspeller.IPA1Analysis/";
	private static String datasetfilename="allFamilies.txt";
	private static Map<String,GeneFamily> families;

	
	public static String getIPA1Directory() {
		return IPA1Directory;
	}
	
	public static void setIPA1Directory(String dir) {
		IPA1Directory=dir;
	}

	public static String getDatasetPath() {
		return IPA1Directory+datasetfilename;
	}

	public static void setDatasetFromFile() throws IOException {
		if (families != null){
			return;
		}
		
		families = new HashMap<String,GeneFamily>();
		BufferedReader in = new BufferedReader(new FileReader(IPA1Directory+datasetfilename));

		GeneFamily gf = new GeneFamily(in);
		while (gf.isInitialized()){
			if (gf.getGenes().size()<=maxGenesInFamily){
				families.put(gf.getFamilyName(),gf);
			}
			gf = new GeneFamily(in);
		}
		in.close();
	}

	 static Map<String, GeneFamily> getFamilies() {
		return families;
	}

	 public static String getActualBS(Suffix s, Sequence seq) {
			if (s.getSequenceID()==0){
				return seq.toString().substring(s.getSequencePosition(),s.getSequencePosition()+maxPatternLength);
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
		
		private static String generateFilename(String theme){
			return IPA1Directory + theme + ".txt";
		}

		public static String getFilenameExperimentalMatches() {
			return generateFilename("IPA1Exp");
		}
		
		public static String getDeNovoFilename(String type,int bls) {
			return generateFilename("DeNovoOCR"+type+bls);
			//return generateFilename("DeNovo"+type+bls);
			
		}
		
		public static String getPMFilename(String type, String motif, int bls) {
			
			return generateFilename("PMsOCR/PMs"+type+motif+bls);
			//return generateFilename("PMs/PMs"+type+motif+bls);
		}
		
		public static String getFilenameExperimentalMatchesWithOrthologs() {
			return generateFilename("IPA1ExpO");
		}
		
		public static String getFilenamePWMMatches() {
			return generateFilename("IPA1PWM");
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
				
				Gene queryGene = new Gene(scanner.next(),species);
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
				
				Gene queryGene = new Gene(scanner.next(),species);
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
		
		
		public static Set<Gene> calculateIntersection(Set<Gene> s1, Set<Gene> s2) {
			Set<Gene> intersection = new HashSet<Gene>();
			for (Gene g : s1){
				if (s2.contains(g)){
					intersection.add(g);
				}
			}
			return intersection;
			
		}

		public static Set<Gene> calculateDifference(Set<Gene> s1, Set<Gene> s2) {
			Set<Gene> diff = new HashSet<Gene>();
			for (Gene g : s1){
				if (!s2.contains(g)){
					diff.add(g);
				}
			}
			for (Gene g : s2){
				if (!s1.contains(g)){
					diff.add(g);
				}
			}
			return diff;
			
		}
}
