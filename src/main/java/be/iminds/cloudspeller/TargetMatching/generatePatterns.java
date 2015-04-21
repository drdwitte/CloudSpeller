package be.iminds.cloudspeller.TargetMatching;

import java.util.HashSet;
import java.util.Set;

public class generatePatterns {

	//used to generate more degenerate patterns starting from the most interesting KN1 vars:
	
	/*
	 	TTTTTT
		GGGGGG
		AAAAAA
		
		TYTYYT	-> 3T3Y 
		NNNNTN	-> 5N1T
		
		GGGGGG
		AAAAAA
		
		YTYTYT	-> 3T3Y
		KWWKTK	-> 3K2W1T
		
		GGGGGG
		AAAAAA
		
		TTTTYY	-> 4T2Y
	*/
			
	
	
	public static void main (String [] args){
		 
		
		Set<String> before = new HashSet<String>();
		Set<String> after = new HashSet<String>();
		
		before.add("TGA");
		
		
		String [] chars = {"TY","TN","G","A","TY","TWK","G","A","TY"};
		
		for (int i=0; i<chars.length; i++){
		
			addChars(before,after,chars[i]);
			reset(before,after);
		}
		
		System.out.println(before.size());
		
		for (String s : before){
			System.out.println(s+"\t15");
		}
		
		
	}
	
	

	private static void addChars(Set<String> before, Set<String> after,
			String chars) {
		
		for (int i=0; i<chars.length(); i++){
			for (String s : before){
				after.add(s+chars.charAt(i));
			}
		}
		
		
	}



	private static void reset(Set<String> before, Set<String> after) {
		before.clear();
		before.addAll(after);
		after.clear();
		
		
	}
	
	
}
