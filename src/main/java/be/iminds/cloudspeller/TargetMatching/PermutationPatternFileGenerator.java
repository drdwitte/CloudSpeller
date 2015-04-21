package be.iminds.cloudspeller.TargetMatching;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;


import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifpermutationgroups.IUPACContent;
import be.iminds.cloudspeller.motifpermutationgroups.MotifContent;

public class PermutationPatternFileGenerator {

	private static final int numberOfPerms=1000;
	/**
	 * Run with java -jar TMGenerator.jar
	 * @param args
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public static void main (String [] args) throws NumberFormatException, IOException{
				
		String filenameTMPatterns = args[0];
		
		BufferedReader in = new BufferedReader(new FileReader(filenameTMPatterns));
		
		String line;
		while ( (line = in.readLine()) !=null)
		{
			String [] splits = line.split("\t");
			String motif = splits[0];
			int minBLS = Integer.parseInt(splits[1]);
			generateTMFile(motif,minBLS);
		}
	}

	public static void generateTMFile(String motif, int BLS) throws IOException {
		
		String output = motif+BLS+".txt";
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		
		
		List<Character> t= new ArrayList<Character>();
		for (int i=0; i<motif.length(); i++){
			t.add(motif.charAt(i));
		}
		Collections.sort(t);
		StringBuilder s = new StringBuilder();
		for (int i=0; i<t.size(); i++){
			s.append(t.get(i));
		}

		MotifContent content = new IUPACContent(s.toString());
		Set<Motif> permutations = content.createPermutationGroup(numberOfPerms);

		permutations.add(new IUPACMotif(motif)); //add query motif!
		
		for (Motif m : permutations){
			writer.write(m.toString()+"\t"+BLS+"\n");
		}
		
		writer.close();
	}		
}
