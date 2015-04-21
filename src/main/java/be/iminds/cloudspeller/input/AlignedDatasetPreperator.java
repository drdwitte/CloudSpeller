package be.iminds.cloudspeller.input;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import be.iminds.cloudspeller.KN1Analysis.KN1Toolbox;

import be.iminds.cloudspeller.toolbox.LineIterator;



public class AlignedDatasetPreperator extends DatasetPreparator {

	public static void main (String [] args) throws IOException {


		int maxGenes=15;
		KN1Toolbox.setDatasetFromFile("allFamilies.txt");
		
		Set<String> fams = KN1Toolbox.getFamilies().keySet();
		System.out.println(fams.size());
		
		
		Set<String> famsProcessed = new HashSet<String>();
		
		String folderWithAlignments = "/home/ddewitte/Desktop/FDRFULLREALANDRANDOM/ShuffledSeqsAlignedAll/";
		String outputDir = "/home/ddewitte/Desktop/FDRFULLREALANDRANDOM/AlignedFams2kbShuffled1/";
		
		//int numberOfGenesFilter = 10;
		
		
		File folder = new File(folderWithAlignments);
		File [] gfFiles = folder.listFiles();
		
		
		Set<GeneFamily> geneFamilies = new HashSet<GeneFamily>();
		int c=0;
		for (File file : gfFiles){
			c++; if (c%1000==0){ System.out.print("*");}
			GeneFamily family =generateFamilyFromABFile(file); 
			
			if (family.getNumberOfGenes()>maxGenes){
				;
			} else {
				geneFamilies.add(family);
				famsProcessed.add(family.getFamilyName());
			}
		}
		
		//adding newicks
		GeneFamily.setGeneralNewick(readNewickFile("speciesTree.txt"));		
		
		for (GeneFamily family : geneFamilies){
			family.generateNewick(paralogBranchLength);
		}
		
		System.out.println("start printing");
		File outputFolder = new File(outputDir);
		outputFolder.mkdir();  
		printGeneFamiliesInGroups(outputDir, geneFamilies, 1);
		
		System.out.println(fams.size() + "\t" + famsProcessed.size());
		for (String fam : fams){
			
			
			if (!famsProcessed.contains(fam)){
				System.out.println("Fam missing: \t"+fam);
			} 
			
		}
	}

	private static GeneFamily generateFamilyFromABFile(File file) throws IllegalArgumentException, IOException {
				
		String familyID = file.getName().split("_")[0];
		//System.out.println(familyID);
		GeneFamily family = new GeneFamily(familyID);
		
		LineIterator iterator = new LineIterator(file.getAbsolutePath());
		
		Gene g = null;
		BaseSequence seq = null;
		while (iterator.hasNext()){
			
			String temp = iterator.next();
			if (temp.length()==0){
				continue;
			}
			
			if (temp.startsWith(">")){
								
				if (g!=null && seq!=null){
					
					family.addGeneSeq(g, seq);
				}
				g = new Gene(temp.substring(1),temp.substring(1,3));
				seq = new BaseSequence("");
			
			} else {
				seq.append(temp);
			}
		}
		
		
		if (g!=null && seq!=null){
			family.addGeneSeq(g, seq);
		}
		
		
	
		return family;
	}
	
	
}
