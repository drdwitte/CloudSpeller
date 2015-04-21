package be.iminds.cloudspeller.motifpermutationgroups;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.alphabets.Alphabet;

	//lazy redundance checker
	/**
	 * For every permutationgroupID we generate the complementary groupID, lexigraphically
	 * first one is not redudant (false) the other is redundant (true)
	 */
public class RedundanceChecker { 
		
	private Map<String,Boolean> redundanceMap = new HashMap<String,Boolean>();
	private MotifContentFactory contentFactory;
	private Alphabet alphabet;
	
	public RedundanceChecker(MotifContentFactory contentFactory, Alphabet alphabet){
		this.contentFactory=contentFactory;
		this.alphabet=alphabet;
		
	}

	public boolean isRedundant(MotifContent content){
		String actualKey = contentFactory.createStringRepresentation(content);
		Boolean value = redundanceMap.get(actualKey);
		if (value!=null){
			return value;
		} else {
			
			String complementaryKey = generateComplement(actualKey);
			int comparison = actualKey.compareTo(complementaryKey);
			
			if (comparison < 0){
				redundanceMap.put(actualKey,false);
				redundanceMap.put(complementaryKey,true);
				return false;
			} else if (comparison > 0){
				redundanceMap.put(actualKey,true);
				redundanceMap.put(complementaryKey,false);
				return true;
			} else {//==
				redundanceMap.put(actualKey,false);
				return false;
			}
			
			
		}
	}

	private String generateComplement(String s) {
				
		char[] chars = s.toCharArray();
	    
		for (int i=0; i<chars.length; i++){
			chars[i]=alphabet.getComplement(chars[i]);
		}
		
		Arrays.sort(chars);
	    return new String(chars);
	}
}

