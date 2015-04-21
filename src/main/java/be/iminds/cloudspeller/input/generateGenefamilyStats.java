package be.iminds.cloudspeller.input;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class generateGenefamilyStats {

	
	private static Map<String,String> genePrefixOrgMap = new HashMap<String, String>(); 

	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		//read different sequence files
		String bdiFile="bdi_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS1=DatasetPreparator.readGeneSequenceFile(bdiFile,"BD");
		System.out.println(bdiFile+" contains "+gS1.size()+" sequences");
		genePrefixOrgMap.put("BD","BD");
			
		String sbiFile="sbi_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS2=DatasetPreparator.readGeneSequenceFile(sbiFile,"SB");
		System.out.println(sbiFile+" contains "+gS2.size()+" sequences");
		genePrefixOrgMap.put("SB","SB");	
		
		String osaFile="osa_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS3=DatasetPreparator.readGeneSequenceFile(osaFile,"OS");
		System.out.println(osaFile+" contains "+gS3.size()+" sequences");
		genePrefixOrgMap.put("OS","OS");	
		
		String zmaFile="zma_plaza2_5_2kbupstream.tfa";
		Map<Gene,String> gS4=DatasetPreparator.readGeneSequenceFile(zmaFile, "ZM");
		System.out.println(zmaFile+" contains "+gS4.size()+" sequences");
		genePrefixOrgMap.put("ZM","ZM");
		
		//join gene sequence maps
		Map<Gene,String> allGeneSeqs=new HashMap<Gene,String>();
		allGeneSeqs.putAll(gS1);
		allGeneSeqs.putAll(gS2);
		allGeneSeqs.putAll(gS3);
		allGeneSeqs.putAll(gS4);
		
		
		Map<String,Set<Gene>> geneFams = DatasetPreparator.readSelectionFile("iORTHO_2evid_osa.selection",genePrefixOrgMap);
		System.out.println(geneFams.size());
		System.out.println("Gene family \t #genes \t #chars");
		for (Map.Entry<String,Set<Gene>> e : geneFams.entrySet()){
			
			StringBuilder s = new StringBuilder();
			
			int nChars = 0;
			boolean allSeqsFound = true;
			for (Gene g : e.getValue()){
				String seq = allGeneSeqs.get(g);
				if (seq==null){ //sequence missing
					allSeqsFound=false;
					break;
				} else {
					nChars+=seq.length();
				}
			}
			
			if (allSeqsFound){
			
				s.append(e.getKey());
				s.append("\t");
				s.append(e.getValue().size());
				s.append("\t");
				s.append(nChars);
				
				System.out.println(s.toString());
			}
			
			
			
		
		}
		
		return;
	}

}
