package be.iminds.testcloudspeller;

import static org.junit.Assert.*;
import be.iminds.cloudspeller.indexing.FullMotifMatch;
import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.IndexStructure;
import be.iminds.cloudspeller.indexing.NoDecorationFactory;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.indexing.PatternBLSPair;
import be.iminds.cloudspeller.indexing.SequenceIDSet;
import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.indexing.SuffixInterpreter;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.Sequence;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;

import org.junit.Test;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;
import be.iminds.cloudspeller.phylogenetics.BLSCalculatorFactory;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculatorFactory;

import be.iminds.cloudspeller.KN1Analysis.BolducAnalyzer;
import be.iminds.cloudspeller.KN1Analysis.DeNovoRecord;
import be.iminds.cloudspeller.KN1Analysis.DegMatchesProcessor;
import be.iminds.cloudspeller.KN1Analysis.KN1Toolbox;
import be.iminds.cloudspeller.KN1Analysis.RecF;
import be.iminds.cloudspeller.KN1Analysis.RecS;

public class KN1AnalysisTest {

	//De Novo testing: ConfidenceChars -> patternBLS -> pmOutput -> 
	
	@Test
	/**
	 * Simply run the KN1PredictionsDeNovo on a 1 be.iminds.cloudspeller.output file and check if the
	 * extraction is succesful by random sampling the file for 20 motifs 
	 * the extraction of 20 motifs is done using grep with TGA..GA..GA. regex
	 */
	public void testDeNovoFromConfidenceChartsToPatternBLS() throws IOException {
		
		//read the be.iminds.cloudspeller.output of the algorithm
		Set<PatternBLSPair> pairs = new HashSet<PatternBLSPair>();
		
		String pmFile="unitTestPMfile.txt";
		DegMatchesProcessor processor = new DegMatchesProcessor(null);
		processor.reworkKN1Predictions("outputUnitTestSifter/part-r-00000", "unitTestTemp.txt",KN1Toolbox.generateKN1PWM());
		DegMatchesProcessor.generateInputFileForDBPatternMatcher("unitTestTemp.txt",pmFile);
		
		BufferedReader in = new BufferedReader(new FileReader(pmFile));
	    String patternBLS;
	    while ((patternBLS = in.readLine()) != null) {
	    	pairs.add(new PatternBLSPair(patternBLS));
	    }
	    in.close();
		  
	    /* RANDOM SAMPLE WITH GREP
	    TGACAGAATGAC	15
	    TGAGMGAKCGAW	50
	    TGAWMGAKCGAG	50
	    TGASMGAYCGAT	50
	    TGACTGASYGAM	15
	  	TGACTGAMYGAS	50
	    TGASYGAMTGAC	15
	    TGAMTGASYGAC	50
	    TGAWRGATCGAK	50
	    TGATRGAWCGAK	50
	    TGAGSGAGMGAG	15
	    TGAGGGAKWGAG	50
	    TGAGWGAGKGAG	15
	    TGAWKGAGGGAG	50
	    TGAWGGAGKGAG	50
    	TGAWGGAGGGAK	50
	    TGATKGAYGGAN	15
	    
	    TGATNGAKYGAG	50
	    TGAYNGATGGAK	50
	    TGAKNGATGGAY	50
	     */
	    
	    String [] patterns = {"TGACAGAATGAC","TGAGMGAKCGAW","TGAWMGAKCGAG","TGASMGAYCGAT","TGACTGASYGAM",
	    					  "TGACTGAMYGAS","TGASYGAMTGAC","TGAMTGASYGAC","TGAWRGATCGAK","TGATRGAWCGAK",
	    					  "TGAGSGAGMGAG","TGAGGGAKWGAG","TGAGWGAGKGAG","TGAWKGAGGGAG","TGAWGGAGKGAG",
	    					  "TGAWGGAGGGAK","TGATKGAYGGAN","TGATNGAKYGAG","TGAYNGATGGAK","TGAKNGATGGAY"};
	    
	    int [] blsScores = {15,50,50,50,15,50,15,50,50,50,15,50,15,50,50,50,15,50,50,50};
	    String pattern;
	    int bls;
	    
	    
	    for (int i=0; i<patterns.length; i++){
	       	pattern = patterns[i]; bls = blsScores[i];
	       	assertTrue(pairs.contains(new PatternBLSPair(pattern,bls)));
	    }
	}
	
	@Test
	public void testDeNovoDistributedPatternMatching(){
		//TODO tested manually results tend to agree with confidencegraphs
	}
	
	@Test
	//see also manualInnerJoinFilteringUnitTest.txt
	public void testInnerJoin() throws IOException{
		
		String outputPM = "PMKN1Variants6SeptemberUnitTest.txt";

		Set<DeNovoRecord> kn1VariantsAbovePWMThreshold;
		kn1VariantsAbovePWMThreshold = DegMatchesProcessor.selectDeNovoRecords(-70,"unitTestTemp.txt"); //92
		assertEquals(92, kn1VariantsAbovePWMThreshold.size());
		kn1VariantsAbovePWMThreshold = DegMatchesProcessor.selectDeNovoRecords(-50,"unitTestTemp.txt"); //76
		assertEquals(76, kn1VariantsAbovePWMThreshold.size());
		kn1VariantsAbovePWMThreshold = DegMatchesProcessor.selectDeNovoRecords(-40,"unitTestTemp.txt"); //37
		assertEquals(37, kn1VariantsAbovePWMThreshold.size());
		kn1VariantsAbovePWMThreshold = DegMatchesProcessor.selectDeNovoRecords(-30,"unitTestTemp.txt"); //33
		assertEquals(33, kn1VariantsAbovePWMThreshold.size());
		kn1VariantsAbovePWMThreshold = DegMatchesProcessor.selectDeNovoRecords(-10,"unitTestTemp.txt"); //8
		assertEquals(8, kn1VariantsAbovePWMThreshold.size());
		
		Map<Gene,Set<RecS>> records = DegMatchesProcessor.innerJoinDeNovoAndPMRecords(kn1VariantsAbovePWMThreshold,outputPM);
		
		//note first motif's maize gene occurs in 2 different families!
		//fourth motif has a maize gene which is in 3 genefamilies
		//third last motif has many duplicates (4 twofolds, 1 threefold)
		//second last has 2 2fold duplates
		//last has 1 twofold
		String [] bestMotifs = {"TGACAGAATGAC","TGAKRGATWGAC", "TGATKGAWRGAC",	"TGAYGGATKGAN",
		"TGAYKGATGGAN",	"TGATGGAYKGAN", "TGATKGAYGGAN",	"TGAYGGACKGAA"};
		
		int [] numMaizeGenesManualCheck = {2,0,1,5,7,21,22,2};
		
		
		for (int i=0; i<bestMotifs.length; i++){
			
			int count = 0;
			for (Map.Entry<Gene,Set<RecS>> e : records.entrySet()){
				for (RecS rec : e.getValue()){
					if (rec.getBindingSite().equals(bestMotifs[i])){
							count++;
					}
				}
			}
			assertEquals(numMaizeGenesManualCheck[i],count);
			
		}
	}
	

	
	
	
	//BOLDUC testing: test if SbolducFiltered and SbolducO are correct
	
	/**
	 * Only test 100 first entries in Bolduc
	 */
	@Test
	public void testBolducMaizeGenesWithReference() throws IOException{
		
		KN1Toolbox.setDatasetFromFile("allFamiliesUnitTest.txt");
		BolducAnalyzer bolducAnalyzer = new BolducAnalyzer();
		Set<Gene> bolducGenes = bolducAnalyzer.getBolducGenesFromFile("KN1TargetsAllFromBolducUnitTest.txt");
		Map<Gene,Set<String>> bolducFamilyMap = bolducAnalyzer.getFamilyOfGenes(bolducGenes);
		
		Map<Gene, Set<RecS>> tableSmall = bolducAnalyzer.getBolducRecordsWithReferenceInPromoter(bolducFamilyMap,KN1Toolbox.referenceMotif);
		System.out.println(tableSmall.size());
		Map<String,GeneFamily> familieMap =  KN1Toolbox.getFamilies();
		
		//checking tableSmall
		for (Map.Entry<Gene, Set<RecS>> e : tableSmall.entrySet()){
			Sequence seq = getSequence(familieMap, e.getKey());
			for (RecS record : e.getValue()){
				
				String forward = record.getBindingSite();
				String reverse = KN1Toolbox.getComplementOf(forward);
				boolean cond = seq.toString().contains(forward) || seq.toString().contains(reverse); 
				
				if (!cond){
					System.out.println("Problem with: "+ e.getKey());
					System.out.println(seq);
					System.out.println(record);
				}
				assertTrue(cond);
			}
		}
		
		//checking orthologs
		Map<Gene, Set<RecF>> tableFull  = bolducAnalyzer.getBolducRecordsWithReferenceInPromoterAndOrthologs(bolducFamilyMap,KN1Toolbox.referenceMotif);

		for (Map.Entry<Gene, Set<RecF>> e : tableFull.entrySet()){
			
			for (RecF record : e.getValue()){
				Sequence seq = getSequence(familieMap, new Gene(record.getGene(),record.getGene().substring(0,2)));
				String forward = record.getBindingSite();
				String reverse = KN1Toolbox.getComplementOf(forward);
				boolean cond = seq.toString().contains(forward) || seq.toString().contains(reverse); 
				
				if (!cond){
					System.out.println("Problem with: "+ e.getKey());
					System.out.println(seq);
					System.out.println(record);
				}
				assertTrue(cond);
			}
		}
		
		Set<String> keys = getFamilies(familieMap, new Gene("ZM01G11980","ZM"));
		
		System.out.println(keys);
		
		
		for (String s : keys){
			
			GeneFamily g = familieMap.get(s);
			
			for (int i=0; i<g.getGenes().size(); i++){
				List<Suffix> matches = KN1Toolbox.giveMatchesOfPatternInSequence("TGANNGANNGAN", g.getSequences().get(i));
				int numMatches;
				if (matches == null){
					numMatches = 0;
				} else {
					numMatches = matches.size();
				}
				System.out.println(g.getGenes().get(i)+"\t"+numMatches);
				
			}
		}
		
		
		
	}
	
	private Sequence getSequence(Map<String,GeneFamily> familieMap, Gene gene){
		for (Map.Entry<String,GeneFamily> e : familieMap.entrySet()){
			for (int i=0; i<e.getValue().getNumberOfGenes(); i++){
				if (e.getValue().getGenes().get(i).equals(gene)){
					return e.getValue().getSequences().get(i);
				}
			}
		}
		return null;
	}
	
	private Set<String> getFamilies(Map<String,GeneFamily> familieMap, Gene gene){
		Set<String> familiesWithGene = new HashSet<String>();
		for (Map.Entry<String,GeneFamily> e : familieMap.entrySet()){
			for (int i=0; i<e.getValue().getNumberOfGenes(); i++){
				if (e.getValue().getGenes().get(i).equals(gene)){
					familiesWithGene.add(e.getKey());
				}
			}
		}
		return familiesWithGene;
	}
	
	@Test
	public void testDBPMWithWeirdMatch(){
		int bls = 15;
		GeneFamily gf = KN1Toolbox.getFamilies().get("iORTHO010829");
		
		GSTFactory factory = new GSTFactory(12,true,new NoDecorationFactory());				
		
		IndexStructure iS = factory.createIndexStructure(gf.getSequences());
		ConservationScoreCalculatorFactory calculatorFac = new BLSCalculatorFactory(new BLS(0));
		BLSCalculator calculator = (BLSCalculator) calculatorFac.createCalculator(gf);
		
		SuffixInterpreter interpreter = new SuffixInterpreter(gf);
		
		String pattern = "TGAKGGATGGAG";
		
		List<Suffix> suffixes = iS.matchExactPattern(new IUPACMotif(pattern));
		BLS cutoffScore = new BLS(bls);
				
		if (suffixes==null){
			System.out.println("No Matches");
		}
			
		NodeDecoration nodeInfo = new SequenceIDSet(gf.getNumberOfGenes());
		nodeInfo.processSuffixes(suffixes);
		BLS score = (BLS) calculator.calculateScore(nodeInfo);
		
		
		System.out.println(score);
		System.out.println(gf.getNewick());
		
		
		if (score == null){ //score > BLSMIN
			System.out.println("Score below threshold");;
		}
		if (score.compareTo(cutoffScore)<0){ //score >=patternBLS
			System.out.println("Score below threshold");;
		}
							
		Set<FullMotifMatch> uniqueMatches = new HashSet<FullMotifMatch>();
				
		for (Suffix s : suffixes){
			FullMotifMatch match = interpreter.translateSuffix(pattern,s,bls);
			uniqueMatches.add(match);
		}
			
			
		for (FullMotifMatch m : uniqueMatches){
			
			System.out.println(m);
			
		}
			
			
			
		
	}
	
	
	@Test
	public void testRecordSelectionInBolduc(){
		
	}
	
	
	//Venn Diagram testing
	@Test
	public void testGenerateVennDiagram(){
		
	}
	
	@Test
	public void testSingleMotifOverlap(){
		
	}
}
