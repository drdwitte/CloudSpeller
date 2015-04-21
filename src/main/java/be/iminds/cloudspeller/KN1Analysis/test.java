package be.iminds.cloudspeller.KN1Analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

public class test {

	public static Map<Character,String> degMap = new HashMap<Character,String>();
	
	public static void init(){
		degMap.put('A',"AMWRN");
		degMap.put('C',"CYSMN");
		degMap.put('G',"GRSKN");
		degMap.put('T',"TYWKN");
		
		degMap.put('R',"RN");
		degMap.put('Y',"YN");
		degMap.put('S',"SN");
		degMap.put('W',"WN");
		degMap.put('K',"KN");
		degMap.put('M',"MN");
		
		degMap.put('N',"N");
	}
	
	public static void main (String [] args){

		init();
		String regex = "TGA..GA..GA.";
		String motif = "TGATYGATKGAT" ;
		List<StringBuilder> candMotifs = createAllDegMatchesWithRegExp(regex,motif);
		
		
		//PWM pwm = KN1Toolbox.generateKN1PWM();

		SortedSet<MotifScore> motifScores = new TreeSet<MotifScore>();
		DegMatchesProcessor processor = new DegMatchesProcessor(null);
		processor.initializeBSMap(regex,KN1Toolbox.generateKN1PWM());
		for (int i=0; i<candMotifs.size(); i++){
			String candMotif = candMotifs.get(i).toString();
			double score = processor.calculatescoreWorstBS(candMotif);
			motifScores.add(new MotifScore(candMotif, score));
		}
		
		for (MotifScore m : motifScores){
			System.out.println(m);
		}
		
	}
	
public static List<StringBuilder> createAllDegMatchesWithRegExp(String regex, String str) {
		
		List<StringBuilder> oldA = new ArrayList<StringBuilder>();
		List<StringBuilder> newA = new ArrayList<StringBuilder>();
		
		oldA.add(new StringBuilder(""));
		
		for (int i=0; i<regex.length(); i++){
			char c = regex.charAt(i);
			
			if (c == '.'){
				for (StringBuilder b : oldA){
					
					String alph = degMap.get(str.charAt(i));
					
					for (int j=0; j<alph.length(); j++){
						
						StringBuilder copy = new StringBuilder(b);
						copy.append(alph.charAt(j));
						newA.add(copy);
					}
				}
				oldA.clear();
				oldA.addAll(newA);
				
				newA.clear();
				//oldA = newA;
				
			} else {
				for (StringBuilder b : oldA){
					b.append(c);
				}
			}
		}
		
		return oldA;
	}

	
	
}
