package be.iminds.cloudspeller.driver;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import be.iminds.cloudspeller.indexing.BitSetDecorationFactory;
import be.iminds.cloudspeller.indexing.GSTFactory;
import be.iminds.cloudspeller.indexing.IndexStructureFactory;
import be.iminds.cloudspeller.indexing.NodeDecorationFactory;

import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.motifalgorithms.*;
import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.phylogenetics.BLSCalculator;
import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

public class RetrieveStatsForASingleGF {

	private static final String baseChars="ACGT";

	private static final int kmin=6;
	private static final int maxDepth = 12;
	private static MotifExtractor extractor;
	private static IUPACType iupacType = IUPACType.FULL;

	private static List<GeneFamily> gfs = new ArrayList<>();

	public static void initializeSearchAlgorithm(DiscoveryAlgorithm alg, GeneFamily gf, int kmax, int maxDegPos){
		//int kmax=12;

		//a)

		boolean withReverseComplements=true;
		NodeDecorationFactory nodeDecoFac= new BitSetDecorationFactory();
		IndexStructureFactory iFac = new GSTFactory(maxDepth, withReverseComplements, nodeDecoFac);
		alg.setDataStructure(iFac);

		//b) 


		Alphabet alph = new IUPACAlphabet(iupacType);
		MotifSearchSpace searchSpace = new MotifSearchSpace(kmin,kmax,maxDegPos,alph);

		alg.setSearchSpace(searchSpace);

		//c)

		BLSCalculator calculator=new BLSCalculator(gf);
		calculator.setCutoff(new BLS(BLS.MIN));
		alg.setConservationScoreCalculator(calculator);

		//d)

		extractor = new DevNullContainer();
		//extractor = new MotifContainer();

		alg.setMotifExtractor(extractor);


	}

	public static void main(String [] args) throws IOException {

		int eMin =3;
		int eMax = 3;

		int lengthMax=kmin;

		String dirPrefix = "/home/ddewitte/Desktop/bioinformaticsPHD/";

		String fileAverageMotifs = dirPrefix + "averageMotifs.tsv";
		String fileAverageDiscoveryTime = dirPrefix + "averageTime.tsv";
		String fileAllSamples = dirPrefix + "samples.tsv";

		BufferedWriter wM = new BufferedWriter(new FileWriter(new File(fileAverageMotifs)));
		BufferedWriter wT = new BufferedWriter(new FileWriter(new File(fileAverageDiscoveryTime)));
		BufferedWriter wS = new BufferedWriter(new FileWriter(new File(fileAllSamples)));


		int [] blsTresholds = {15,50,60,70,90,95};
		BLS.initializeBLSConstants(blsTresholds); 
		
		String dirName = dirPrefix + "performanceTestset/groupOrtho";

		int nGroups = 10;
		int gfsPerFile = 10;

		for (int i=1; i<=nGroups; i++){

			BufferedReader in = new BufferedReader(new FileReader(new File(dirName+i+".txt")));

			for (int j=0; j<gfsPerFile; j++) {
				gfs.add(new GeneFamily(in));

			}
		}

		for (int l=kmin; l<=lengthMax; l++) {

			for (int e=eMin; e<=eMax; e++) {
				List<Integer> nMotifs = new ArrayList<>();
				List<Double> timePerFam = new ArrayList<>();

				for (int i = 0; i < gfs.size(); i++) {
					GeneFamily gf = gfs.get(i);
					//System.out.println(gf);
					DiscoveryAlgorithm alg = new DeNovoExactDiscoveryAlgorithm(gf.getSequences());

					initializeSearchAlgorithm(alg, gf, l, e);

					long start = System.nanoTime();
					alg.runDiscovery(new IUPACFactory(iupacType));
					long stop = System.nanoTime();


					double numMilliSecs = (stop - start) / 1000000;
					int numberOfMotifs = extractor.getNumberOfMotifsExtracted();

					//alg lopen: aantal motieven / tijd voor discovery


					wS.write(i + "\t" + e + "\t" + l + "\t" + numberOfMotifs + "\t" + numMilliSecs); wS.newLine();
					nMotifs.add(numberOfMotifs);
					timePerFam.add(numMilliSecs);
					System.out.println(i + "\t" + e + "\t" + l + "\t" + numberOfMotifs + "\t" + numMilliSecs);
				}

				double avgMotifs = 0;
				double avgTime = 0.0;

				for (int i = 0; i < nMotifs.size(); i++) {
					avgMotifs += nMotifs.get(i) * 1.0 / nMotifs.size();
					avgTime += timePerFam.get(i) * 1.0 / nMotifs.size();
				}

				wM.write(e+"\t"+l+"\t"+avgMotifs); wM.newLine();
				wT.write(e+"\t"+l+"\t"+avgTime); wT.newLine();
			}
		}



		wM.close();
		wS.close();
		wT.close();









	}



}
