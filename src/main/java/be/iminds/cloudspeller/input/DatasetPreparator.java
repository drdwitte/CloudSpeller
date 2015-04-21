package be.iminds.cloudspeller.input;

import org.apache.avro.generic.GenericData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.*;

public class DatasetPreparator {

	protected static final double paralogBranchLength = 0.000001;
	protected static Map<String,String> genePrefixOrgMap = new HashMap<String, String>(); 
	
	public static String readNewickFile(String filename) throws IOException{
		BufferedReader in=new BufferedReader(new FileReader(filename));
		String newick = in.readLine(); 
		in.close();
    	return newick;
	}
	
	/**
	 * Line to process:
	 * 	#		query_gene	label	ortho_group		#genes	#species	species_list		inparalogs		gene_list
		COUNT	OS01G01010	valid	iORTHO000001	4		4			bdi osa sbi zma		OS01G01010		BD1G74660 OS01G01010 SB03G009250 ZM03G07960
	 */
	
	public static Map<String,Set<Gene> > readSelectionFile(String filename, Map<String,String> prefixOrgMap) throws IOException {
		Map<String,Set<Gene> > geneFams= new HashMap<String,Set<Gene> >(); 
		BufferedReader in=new BufferedReader(new FileReader(filename));
		in.readLine(); //skip first line
		String line;
		while ((line = in.readLine()) != null) {
			if (line.length()==0){ //ignore empty lines
				continue;
			}
			
			Scanner scan=new Scanner(line);
			scan.next(); scan.next(); scan.next(); //COUNT, query_gene, lqbel
			String orthoGroup=scan.next();
			scan.nextInt(); //number of genes
			int numberOfSpecies=scan.nextInt();
			
			for (int i=0; i<numberOfSpecies; i++){ //ignore species list
				scan.next();
			}
			
			
			
			Set<Gene> genesInFamily= new HashSet<Gene>();
			
			while (scan.hasNext()){ //process list of genes and inparalogs (contains duplicates)
				String geneID=scan.next();
				String org="";
				//identify organism				
				for (Map.Entry<String,String> entry : prefixOrgMap.entrySet()){
					if (geneID.startsWith(entry.getKey())){
						org=entry.getValue();
						break;
					}
				}
				
				if (org.length()!=0){
					genesInFamily.add(new Gene(geneID,org));
				}
			}
			
			scan.close();
						
			geneFams.put(orthoGroup,genesInFamily);
		}
		
		in.close();
		
		return geneFams;
	}
	
	public static Map<Gene,String> readGeneSequenceFile(String filename, String species) throws IOException {
		Map<Gene,String> geneSeqMap=new HashMap<Gene,String>();
		BufferedReader in;
		try {
			in=new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException f){
			System.out.println("Warning: sequence file not found: "+filename);
			return null;
		}
		String line;
		while ((line=in.readLine())!=null){
			String sequence=in.readLine();
			Gene g= new Gene(line.substring(1),species);
			geneSeqMap.put(g,sequence);
		}
		in.close();
		return geneSeqMap;
	}
	
	
	@SuppressWarnings("unused")
	private static void printGeneFamiliesToFile(String filename, Set<GeneFamily> geneFamilies) throws FileNotFoundException, UnsupportedEncodingException {
		try{
		  
		  FileWriter fstream = new FileWriter(filename);
		  BufferedWriter out = new BufferedWriter(fstream);
		  
		  Iterator<GeneFamily> familyIterator = geneFamilies.iterator();

		  while (familyIterator.hasNext()){
			  out.write(familyIterator.next().toString());
		  }
		  
		  
		  out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
		
	}
	
	@SuppressWarnings("unused")
	private static void printGeneFamiliesToTheirOwnFile(String mapName, Set<GeneFamily> geneFamilies, int max) throws FileNotFoundException, UnsupportedEncodingException {
		
		Iterator<GeneFamily> familyIterator = geneFamilies.iterator();
		String extension = ".txt";
		int numFiles=0;
		while (familyIterator.hasNext() && numFiles < max){
			
			GeneFamily gf = familyIterator.next();
			try{
				
				
				String filename = mapName + gf.getFamilyName() + extension;
		  
				FileWriter fstream = new FileWriter(filename);
				BufferedWriter out = new BufferedWriter(fstream);
				out.write(gf.toString());
				out.close();
				
			}catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}
			
			numFiles++;
		}
	}
	
	
	protected static void printGeneFamiliesInGroups(String mapName,
			Set<GeneFamily> geneFamilies, int numFamPerFile) throws IOException {

		
		Iterator<GeneFamily> familyIterator = geneFamilies.iterator();
		String extension = ".txt";
		
		String prefix = "groupOrtho";
		//filename groupOrtho + size + _ + ID + .txt
		
		String filename = mapName+prefix+"1"+extension;
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		int counter=0;
		while (familyIterator.hasNext()){
			GeneFamily gf = familyIterator.next();
			out.write(gf.toString());
			counter++;
			
			if (counter%numFamPerFile == 0){
				out.close();
				int id = counter/numFamPerFile +1;
				filename = mapName+prefix+id+extension;
				out = new BufferedWriter(new FileWriter(filename));
			}
	
		}
		
		out.close();
	
	}
	
	
	public static Map<Gene,String> readAllSequenceFiles(String dir) throws IOException{
		//read different sequence files
		String bdiFile="bdi_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS1=readGeneSequenceFile(dir+bdiFile,"BD");
		System.out.println(bdiFile+" contains "+gS1.size()+" sequences");
		genePrefixOrgMap.put("BD","BD");
			
		String sbiFile="sbi_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS2=readGeneSequenceFile(dir+sbiFile,"SB");
		System.out.println(sbiFile+" contains "+gS2.size()+" sequences");
		genePrefixOrgMap.put("SB","SB");	
		
		String osaFile="osa_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS3=readGeneSequenceFile(dir+osaFile,"OS");
		System.out.println(osaFile+" contains "+gS3.size()+" sequences");
		genePrefixOrgMap.put("OS","OS");	
		
		String zmaFile="zma_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS4=readGeneSequenceFile(dir+zmaFile, "ZM");
		System.out.println(zmaFile+" contains "+gS4.size()+" sequences");
		genePrefixOrgMap.put("ZM","ZM");	
		
		//join gene sequence maps
		Map<Gene,String> allGeneSeqs=new HashMap<Gene,String>();
		allGeneSeqs.putAll(gS1);
		allGeneSeqs.putAll(gS2);
		allGeneSeqs.putAll(gS3);
		allGeneSeqs.putAll(gS4);
			
		System.out.println("Total Number of geneSeqs: "+allGeneSeqs.size());
		
		return allGeneSeqs;
	}

	public static Map<Gene,String> readAllSequenceFilesRandomized(String dir) throws IOException{
		//read different sequence files
		String bdiFile="bdi_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS1=readGeneSequenceFile(dir+bdiFile,"BD");
		System.out.println(bdiFile+" contains "+gS1.size()+" sequences");
		gS1 = randomizeGSRelation(gS1);
		genePrefixOrgMap.put("BD","BD");

		String sbiFile="sbi_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS2=readGeneSequenceFile(dir+sbiFile,"SB");
		System.out.println(sbiFile+" contains "+gS2.size()+" sequences");
		gS2 = randomizeGSRelation(gS2);

		genePrefixOrgMap.put("SB","SB");

		String osaFile="osa_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS3=readGeneSequenceFile(dir+osaFile,"OS");
		System.out.println(osaFile+" contains "+gS3.size()+" sequences");
		gS3 = randomizeGSRelation(gS3);

		genePrefixOrgMap.put("OS","OS");

		String zmaFile="zma_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS4=readGeneSequenceFile(dir+zmaFile, "ZM");
		System.out.println(zmaFile+" contains "+gS4.size()+" sequences");
		gS4 = randomizeGSRelation(gS4);

		genePrefixOrgMap.put("ZM","ZM");

		//join gene sequence maps
		Map<Gene,String> allGeneSeqs=new HashMap<Gene,String>();
		allGeneSeqs.putAll(gS1);
		allGeneSeqs.putAll(gS2);
		allGeneSeqs.putAll(gS3);
		allGeneSeqs.putAll(gS4);

		System.out.println("Total Number of geneSeqs: "+allGeneSeqs.size());

		return allGeneSeqs;
	}

	private static Map<Gene,String> randomizeGSRelation(Map<Gene, String> gS) {
		System.out.println("randomize");
		List<Gene> original = new ArrayList<>();
		original.addAll(gS.keySet());
		List<Gene> shuff = new ArrayList<>();
		shuff.addAll(gS.keySet());

		Collections.shuffle(shuff);

		Map<Gene, String> gShuffled = new HashMap<>();

		for (int i=0; i<original.size(); i++){
			gShuffled.put(shuff.get(i),gS.get(original.get(i)));
		}

		return gShuffled;

	}

	public static Set<GeneFamily> createGeneFamilies(Map<String,Set<Gene>> geneFams, Map<Gene,String> allGeneSeqs){
		//create set of genefamilies
		Set<GeneFamily> geneFamilies = new HashSet<GeneFamily>();
		int numSeqsMissing=0;
		for (Map.Entry<String, Set<Gene> > entry : geneFams.entrySet()){
			
			GeneFamily gf=new GeneFamily(entry.getKey());
			Set<Gene> genes=entry.getValue();
			boolean allSeqsFound=true;
			
			for (Gene g : genes){
				
				String seq=allGeneSeqs.get(g);
		
				if (seq!=null){
					gf.addGeneSeq(g, new BaseSequence(seq));
				} else {
					numSeqsMissing++;
					allSeqsFound=false;
				}
			}
			
			if (allSeqsFound){ //create gene family if all gene seqs are in the dataset
				gf.generateNewick(paralogBranchLength);
				//System.out.println(gf.getNewick());
				geneFamilies.add(gf);
			}
		}
		
		System.out.println("Number of missing gene seqs: "+numSeqsMissing);
		
		return geneFamilies;
	}
	


	public static void generateRegularDataset() throws IOException{
		String dir = "/home/ddewitte/Desktop/bioinformaticsPHD/OrigFormaatDataset/";
		Map<Gene,String> allGeneSeqs = readAllSequenceFiles(dir);

		Map<String,Set<Gene>> geneFams = readSelectionFile(dir+"iORTHO_2evid_osa.selection",genePrefixOrgMap);

		System.out.println("Number of gene families in selection file: "+geneFams.size());


		GeneFamily.setGeneralNewick(readNewickFile(dir+"speciesTree.txt"));





		Set<GeneFamily> geneFamilies = createGeneFamilies(geneFams, allGeneSeqs);





		System.out.println("Printing to file...");
		//printGeneFamiliesToTheirOwnFile("MonocotFamilies2kb_all/", geneFamilies, 18000);

		int numFamPerFile=10;
		printGeneFamiliesInGroups(dir+"MonocotFamiliesAgg10PBL1em6/",geneFamilies,numFamPerFile);
		System.out.println("Done! Bye!");



	}

	public static void generateRandomizedGSRelationshipDataset() throws IOException {
		String dir = "/home/ddewitte/Desktop/bioinformaticsPHD/OrigFormaatDataset/";

		//RANDOMIZED RELATION!!
		Map<Gene,String> allGeneSeqs = readAllSequenceFilesRandomized(dir);


		Map<String,Set<Gene>> geneFams = readSelectionFile(dir+"iORTHO_2evid_osa.selection",genePrefixOrgMap);

		System.out.println("Number of gene families in selection file: "+geneFams.size());


		GeneFamily.setGeneralNewick(readNewickFile(dir+"speciesTree.txt"));


		Set<GeneFamily> geneFamilies = createGeneFamilies(geneFams, allGeneSeqs);


		System.out.println("Printing to file...");
		//printGeneFamiliesToTheirOwnFile("MonocotFamilies2kb_all/", geneFamilies, 18000);

		int numFamPerFile=10;
		printGeneFamiliesInGroups(dir+"MonocotFamiliesAgg10PBL1em6/",geneFamilies,numFamPerFile);
		System.out.println("Done! Bye!");



	}






	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		generateRandomizedGSRelationshipDataset();
	}




	



}

