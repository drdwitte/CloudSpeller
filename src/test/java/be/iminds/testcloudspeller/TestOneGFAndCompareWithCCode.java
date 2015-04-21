package be.iminds.testcloudspeller;

import java.io.BufferedWriter;
import java.io.FileWriter;


import be.iminds.cloudspeller.motifalgorithms.DeNovoExactDiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.DevNullContainer;
import be.iminds.cloudspeller.motifalgorithms.DiscoveryAlgorithm;
import be.iminds.cloudspeller.motifalgorithms.MotifSearchSpace;

import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.Motif;

import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

import be.iminds.cloudspeller.indexing.BitSetDecorationFactory;
import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.NodeDecorationFactory;
import be.iminds.cloudspeller.input.BaseSequence;

import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;

import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;
import be.iminds.cloudspeller.phylogenetics.ConservationScore;

public class TestOneGFAndCompareWithCCode {

	public static final String baseChars="ACGT";
	public static final boolean withReverseComplements = true;

	private static GeneFamily family;
	private static NodeDecorationFactory nodeDecoFac;
	private static String newick="((BD1G74660:0.2688,OS01G01010:0.2688):0.0538,(SB03G009250:0.086,ZM03G07960:0.086):0.2366);";

	private static BaseSequence seqBD = new BaseSequence("TACATTCATATTTAGTCAAATTTGAGACAATTAATATGGATCGGAGGAGAAATTCTCAATCCTCTGAATTTTCTTCTTCATGTATTTCTTAAGTCCTGCAGTCCAAACAGTTCCTAAATTTTTGTTGCCACTCCTGTGAAGTCAAGAGCAGCAACTACCATCCAAGTACCGTTCAGTCAGTCGCCTCTCCTGCTCTCCTTCGCAGCAAAACCGGCCCCGGTCCATTTCAACGCATCTCTTCTAGTCAAACCGCTGGACACAACTCCAGCTTTTCGAAATCCCGGAAATACCCCCGCGGCAACCAACAAGCCGTCCTCGTCTCCGCTCTCCCGGCGCTGCCTGCTCTCGACGGCGGCGATCCCCTGCATCTCAGCGGCAGCGATCCACTCCTCTTCCAACGTCACCCCGATATTCCGCCGTGCTCGCAGCTTCCACAAATTCTGGCCGTGTTTGTGCCCACCCGGCACACCATCTAGATCCGCCCCTGAAGTCTCCGATCC");
	private static BaseSequence seqOS = new BaseSequence("CAGATAAACCACACCCACAGGCACCACCGTCCTTGTTGGTAATGAAGAAGACGAGACGACGACTTCCCCACTAGGAAACACGACGGAGGCGGAGATGATCGACGGCGGAGAGAGCTACAGAAACATCGATGCCTCCTGTCCAATCCCCCCATCCCATTCGGTAGTTGGATTGAAGACTACCGAATAAGAGAAGCAGGCAGGCAGACAAACCCTTGAACCAAGGAGTCCTCGCTGAGGAAGCTTTGGATCCACGACGCAGCTATGGCCTCCCCGCCCACCAGGCCGCCAGCCACAACCAGCTGACTAGGTAGGCTTCCTAGGTAGGGATCCCATCCCTTCGATTCCCTACTCCCTCCCCCGATTGATTTGATTTGATTTGATTTAATTCGATTGCCTGCTTTTCAGGTCGCATGCATCATCAGATTTCAATCTCCCTTCGTTCCCTGTCCCTAATCCAATACCAATAGGGAGCAATCAGCTGCTCCTCGACGGCGAGGGAG");
	private static BaseSequence seqSB = new BaseSequence("CTATTAGTACTATACTACGCCGCTAATAATCAGTCGTGCGGCATCGTCGGCGACCTTGCCGCCGGCGACGCGGTTACACAACTAACCCGGGGCTTTTGAGCTACAGCCTTCTCTCGTGGCTTGTGGGGGAGGGCATATCCCCTTCCTCCCTCCCTCTTCGTCGTCGCTGCCCTCTGCTGCGACTGGACTGGCGTGGAGTGGAGATCCATCCACCCGCACTTGCCACGATCTCTTCCTTCCCGTGCCGGTAACCGCTCCCCGAGGTTACTTGGACGGATTGCTCTGGGGGACTGGCTTTCTTTCTGGCAGCCTAATTTAAATCTTCCTCTCTCACCATTCAGTTGGTCGCGGATTTGCTTGCTCTGGTGCGTGCGCCGCCCGGCAACCCAATCCTGGTTGGGTTCGATTTGCTTCGCCTTCGTTCGTCTTCAGCGCTTGAGATCAGTTTAACCGGCGCGGTTTGGTTGTCCGAAGGCTTCTAGAGTCCAAGTGGGACCCCC");
	private static BaseSequence seqZM = new BaseSequence("CTTGCCCCTCCCCGCGATTCTATTCTACGCCGCTATAATCAGTCGTGCGGCGTCGTCGGCGACCTTGCCGCCTGCGGCGCGGACCACGCGGTTATTATACAGCTAACCGGGGCTTTGAGCTACAGCCTTCTCTTGTGGGGAGGCATCCCCGTCCTCTACTTCGTCGCTGCCCTCTGCTGCGCCTGGACTGGCGTGGAGTGGAGATCCATCCACCCGCACTTGCCACGATCTCTTCCTTCTCGTGTCGGTAACCGCCCCTCGAGGTTGCTTGGACGGATTGCACTAGGGACTGGCTTTCTCGCAACCTAATCTAAATCTTCTTCTCTCACCATTCAGAGTTGGTTGCGGATTTGCTTACTTTGGTGCGCCGCCCGGCAATCCAAACAATCCTGCTTGGGTTCGATTCGCTTCGCCTTCGTCTTCAGCGCTCAAGATCAGTTTAGCCGAGGTGCCGGGCGCGGTTTGGCTGTCCGGAGGCTTTTAGAGTTACTGGAACCCCG");
	
	
	public static void main(String [] args) {
		
		//Vooreerst worden BLS.MIN en aantal BLS intervallen (6) vastgelegd in de klasse BLS, 
		//deze klasse is o.a. verantwoordelijk voor de constructie van frequentievectoren
		int [] thresholds = {10,50,60,70,90,95};
		BLS.initializeBLSConstants(thresholds);
		//parameters zoekruimte
		int kmin=6; System.out.println("kmin= "+kmin);
		int kmax=12; System.out.println("kmax= "+kmax);
		
		int degPos1=3; System.out.println("degPos1= "+degPos1);
		int degPos2=3; System.out.println("degPos2= "+degPos2);
		
		
		//Node decoration factory genereert de informatie in de interne knopen van
		//de suffixboom, we kiezen voor de be.iminds.cloudspeller.indexing.BitSetDecoration (een uint_32)
		nodeDecoFac= new BitSetDecorationFactory();
		
		//Er wordt een genfamilie aangemaakt bestaande uit 4 genen + sequenties, 
		//een newick string en een ID
		family=new GeneFamily("iORTHO000001");
		family.setNewick(newick);
			
		family.addGeneSeq(new Gene("BD1G74660","BD"), seqBD);
		family.addGeneSeq(new Gene("OS01G01010","OS"), seqOS);
		family.addGeneSeq(new Gene("SB03G009250","SB"), seqSB);
		family.addGeneSeq(new Gene("ZM03G07960","ZM"), seqZM);
		
		//Eerste testrun met DC alphabet = ACGT + N 
		System.out.println("DC ALPHABET");

		//discovery wordt getest met een zoekruimte woordlengte = [kmin,kmax], aantal ontaarde posities, alphabet
		double timeDC = testDiscovery(new MotifSearchSpace(kmin,kmax,degPos1,new IUPACAlphabet(IUPACType.DONTCARES)));
		
		System.out.println("Tijd: "+timeDC+" ms");
		
		//Tweede testrun met 11-letterig alfabet = ACGT + N + 6 2voudig ontaarde chars
		System.out.println("TWOFOLDSANDN ALPHABET");

		double time2F = testDiscovery(new MotifSearchSpace(kmin,kmax,degPos2,new IUPACAlphabet(IUPACType.TWOFOLDSANDN)));

		System.out.println("Tijd: "+time2F+" ms");
		
		
		
		
		
	}
	
public static double testDiscovery(MotifSearchSpace searchSpace) {
		//hier worden alle factories geintialiseerd, normaal gezien gebeurt dit in testcloudspeller
		//op basis van de settings (nu nog in be.iminds.cloudspeller.driver.frameWorkTest.run()
	
		//er wordt een exact algoritme aangemaakt (exact = geen errors, wel ontaarde karakters)
		DiscoveryAlgorithm alg = new DeNovoExactDiscoveryAlgorithm(family.getSequences());
		
		//BLSCalculator berekent de BLS scores mbv lookup
		BLSCalculator calculator=new BLSCalculator(family);
		//BLS cutoff wordt geset (motieven met BLS < cutoff krijgen een score 'null'
		calculator.setCutoff(new BLS(BLS.MIN));
		alg.setConservationScoreCalculator(calculator);
		
		//maximale diepte van de suffix tree
		int maxDepth=searchSpace.getMaxLength();
		
		//de GST factory is verantwoordelijk voor de constructie van de GST (suffixtree)
		alg.setDataStructure(new GSTFactory(maxDepth, withReverseComplements, nodeDecoFac));
		
		//zoekruimte wordt geset
		alg.setSearchSpace(searchSpace);
		long start=System.nanoTime();
		
		//motif extractor vangt de resultaten op het discovery algoritme, er zijn een aantal
		//varianten: DevNull telt enkel, InstantEmitter wordt gebruikt in de mapper, 
		//MotifContainer vangt alles op in een hashmap
		DevNullContainer extractor = new DevNullContainer();
		alg.setMotifExtractor(extractor);
		
		//de discovery start op en er wordt een ISMonkey gealloceerd die fungeert als een soort
		//iterator structuur voor een willekeurige indexstructuur.
		//Op dit moment wordt de implementatie van ExactGSMonkey gebruikt (dit is een inner class
		//in index.GeneralizedSuffixTree)
		alg.runDiscovery(new IUPACFactory(IUPACType.FULL));
		long stop=System.nanoTime();
		
		//voor details over het discovery algoritme zie be.iminds.cloudspeller.motifalgorithms.DeNovoExactDiscoveryAlgorithm
		//voor de indexstructuur zie be.iminds.cloudspeller.indexing.GeneralizedSuffixTree
		
		System.out.println(((DeNovoExactDiscoveryAlgorithm)alg).getTempNumberOfMotifs()+ " motifs found");
		int nCppMotifs=0;
		
			
		//dit werkt enkel met MotifContainer, hier worden alle motieven die starten of 
		//eindigen met een een N verwijderd gezien dat in de code (voorlopig) niet gebeurt
		
		/*SortedSet<Motif> motifSet = new TreeSet<Motif>();
		
		for (Map.Entry<Motif,ConservationScore> e : extractor.getMotifMap().entrySet()){
			Motif m = e.getKey();
			if (m.toString().startsWith("N")){
				continue;
			}
			if (m.toString().endsWith("N")){
				continue;
			}
			
			motifSet.add(m);
			
			
			nCppMotifs++;
		}*/
		
		System.out.println("Number of c++ motifs found: "+nCppMotifs);
		
		/*String filename = "SearchSpace"+searchSpace.getAlphabet().getAllChars()+".txt";
		Output out = new Output(filename);
		
		
		
		for (Motif m : motifSet){
			out.writeMotif(m,extractor.getMotifMap().get(m));
		}
		
		out.close();*/
		
		return (stop-start)/1000000;
	
	}


	static class Output {
		
		FileWriter fstream;
		BufferedWriter out;
		int counter=0;
		
		public Output(String filename){
			
			// Create file 
			try{
				fstream = new FileWriter(filename);
				out = new BufferedWriter(fstream);
			
			}catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}
		}
		
		public void close(){
			try{
				out.close();
				fstream.close();
			}catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}
		}
		
		public void writeMotif(Motif m, ConservationScore score){
			try{
				out.write((++counter)+ " "+m.toString()+"\t"+score.toString() + "\n");
			}catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}	
		}
	}

	
}

	


//iORTHO000001.txt
//>iORTHO000001
//((BD1G74660:0.2688,OS01G01010:0.2688):0.0538,(SB03G009250:0.086,ZM03G07960:0.086):0.2366);
//4
//BD1G74660	BD
//TACATTCATATTTAGTCAAATTTGAGACAATTAATATGGATCGGAGGAGAAATTCTCAATCCTCTGAATTTTCTTCTTCATGTATTTCTTAAGTCCTGCAGTCCAAACAGTTCCTAAATTTTTGTTGCCACTCCTGTGAAGTCAAGAGCAGCAACTACCATCCAAGTACCGTTCAGTCAGTCGCCTCTCCTGCTCTCCTTCGCAGCAAAACCGGCCCCGGTCCATTTCAACGCATCTCTTCTAGTCAAACCGCTGGACACAACTCCAGCTTTTCGAAATCCCGGAAATACCCCCGCGGCAACCAACAAGCCGTCCTCGTCTCCGCTCTCCCGGCGCTGCCTGCTCTCGACGGCGGCGATCCCCTGCATCTCAGCGGCAGCGATCCACTCCTCTTCCAACGTCACCCCGATATTCCGCCGTGCTCGCAGCTTCCACAAATTCTGGCCGTGTTTGTGCCCACCCGGCACACCATCTAGATCCGCCCCTGAAGTCTCCGATCC
//OS01G01010	OS
//CAGATAAACCACACCCACAGGCACCACCGTCCTTGTTGGTAATGAAGAAGACGAGACGACGACTTCCCCACTAGGAAACACGACGGAGGCGGAGATGATCGACGGCGGAGAGAGCTACAGAAACATCGATGCCTCCTGTCCAATCCCCCCATCCCATTCGGTAGTTGGATTGAAGACTACCGAATAAGAGAAGCAGGCAGGCAGACAAACCCTTGAACCAAGGAGTCCTCGCTGAGGAAGCTTTGGATCCACGACGCAGCTATGGCCTCCCCGCCCACCAGGCCGCCAGCCACAACCAGCTGACTAGGTAGGCTTCCTAGGTAGGGATCCCATCCCTTCGATTCCCTACTCCCTCCCCCGATTGATTTGATTTGATTTGATTTAATTCGATTGCCTGCTTTTCAGGTCGCATGCATCATCAGATTTCAATCTCCCTTCGTTCCCTGTCCCTAATCCAATACCAATAGGGAGCAATCAGCTGCTCCTCGACGGCGAGGGAG
//SB03G009250	SB
//CTATTAGTACTATACTACGCCGCTAATAATCAGTCGTGCGGCATCGTCGGCGACCTTGCCGCCGGCGACGCGGTTACACAACTAACCCGGGGCTTTTGAGCTACAGCCTTCTCTCGTGGCTTGTGGGGGAGGGCATATCCCCTTCCTCCCTCCCTCTTCGTCGTCGCTGCCCTCTGCTGCGACTGGACTGGCGTGGAGTGGAGATCCATCCACCCGCACTTGCCACGATCTCTTCCTTCCCGTGCCGGTAACCGCTCCCCGAGGTTACTTGGACGGATTGCTCTGGGGGACTGGCTTTCTTTCTGGCAGCCTAATTTAAATCTTCCTCTCTCACCATTCAGTTGGTCGCGGATTTGCTTGCTCTGGTGCGTGCGCCGCCCGGCAACCCAATCCTGGTTGGGTTCGATTTGCTTCGCCTTCGTTCGTCTTCAGCGCTTGAGATCAGTTTAACCGGCGCGGTTTGGTTGTCCGAAGGCTTCTAGAGTCCAAGTGGGACCCCC
//ZM03G07960	ZM
//CTTGCCCCTCCCCGCGATTCTATTCTACGCCGCTATAATCAGTCGTGCGGCGTCGTCGGCGACCTTGCCGCCTGCGGCGCGGACCACGCGGTTATTATACAGCTAACCGGGGCTTTGAGCTACAGCCTTCTCTTGTGGGGAGGCATCCCCGTCCTCTACTTCGTCGCTGCCCTCTGCTGCGCCTGGACTGGCGTGGAGTGGAGATCCATCCACCCGCACTTGCCACGATCTCTTCCTTCTCGTGTCGGTAACCGCCCCTCGAGGTTGCTTGGACGGATTGCACTAGGGACTGGCTTTCTCGCAACCTAATCTAAATCTTCTTCTCTCACCATTCAGAGTTGGTTGCGGATTTGCTTACTTTGGTGCGCCGCCCGGCAATCCAAACAATCCTGCTTGGGTTCGATTCGCTTCGCCTTCGTCTTCAGCGCTCAAGATCAGTTTAGCCGAGGTGCCGGGCGCGGTTTGGCTGTCCGGAGGCTTTTAGAGTTACTGGAACCCCG


