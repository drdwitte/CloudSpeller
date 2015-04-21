package be.iminds.cloudspeller.KN1Analysis;

import be.iminds.cloudspeller.indexing.Suffix;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;
import be.iminds.cloudspeller.input.Sequence;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class FindMisalignedMotif {
	private static int motifID=0;
	private static BolducAnalyzer b = new BolducAnalyzer();
	static {
		try {
			KN1Toolbox.setDatasetFromFile("allFamilies.txt");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void processMotif(String m, Set<Gene> targetGenes){
		System.out.println("Top motif "+(++motifID)+": "+m);
		Map<Gene,Set<String>> bolducFamilyMap = b.getFamilyOfGenes(targetGenes);
		Map<String, GeneFamily> fams = KN1Toolbox.getFamilies();


		for (Map.Entry<Gene,Set<String>> e : bolducFamilyMap.entrySet()){
			for (String s : e.getValue()){
				System.out.println(e.getKey()+"\t"+s);
				GeneFamily gf = fams.get(s);
				List<Gene> gs = gf.getGenes();
			
				for (Gene g : gs){
					Sequence seq = gf.getSequence(g);
					List<Suffix> suffs = KN1Toolbox.giveMatchesOfPatternInSequence(m,seq);
					if (suffs==null) continue;
					for (Suffix su : suffs){
						System.out.print(KN1Toolbox.getActualBS(su,seq)+" "+g.getID()+" "+su+",");
					}
					System.out.println();
				}	
			}
			System.out.println("");
		}
	}
	
	public static void main (String [] args) throws IOException{

		String m = "RTCMATCNATCA";
		Set<Gene> targetGenes = new HashSet<Gene>();
		targetGenes.add(new Gene("ZM01G31480","ZM"));
		targetGenes.add(new Gene("ZM07G10870","ZM"));
		targetGenes.add(new Gene("ZM01G13830","ZM"));
		targetGenes.add(new Gene("ZM10G20990","ZM"));
		targetGenes.add(new Gene("ZM01G54560","ZM"));
		targetGenes.add(new Gene("ZM03G28480","ZM"));
		targetGenes.add(new Gene("ZM01G05520","ZM"));
		targetGenes.add(new Gene("ZM01G40560","ZM"));
		
		
		processMotif(m,targetGenes);
		
		targetGenes.clear();
		m="ATCMRTCNATCA";
		targetGenes.add(new Gene("ZM07G10870","ZM"));
		targetGenes.add(new Gene("ZM01G31480","ZM"));
		targetGenes.add(new Gene("ZM01G13830","ZM"));
		targetGenes.add(new Gene("ZM01G54560","ZM"));
		targetGenes.add(new Gene("ZM03G28480","ZM"));
		targetGenes.add(new Gene("ZM01G40560","ZM"));
		
		processMotif(m,targetGenes);
		
		targetGenes.clear();
		m = "TGAYNGATKGAT";
		
		targetGenes.add(new Gene("ZM01G31480","ZM")); 
		targetGenes.add(new Gene("ZM07G10870","ZM"));
		targetGenes.add(new Gene("ZM01G13830","ZM"));
		targetGenes.add(new Gene("ZM01G54560","ZM"));
		targetGenes.add(new Gene("ZM03G28480","ZM"));
		targetGenes.add(new Gene("ZM03G08480","ZM"));
		targetGenes.add(new Gene("ZM01G40560","ZM"));
		
		processMotif(m,targetGenes);
		
	}
}
