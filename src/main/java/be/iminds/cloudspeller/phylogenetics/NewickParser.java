package be.iminds.cloudspeller.phylogenetics;

import java.util.ArrayList;
import java.util.List;

public class NewickParser
{

	
	/**
	 * "(str)" -> "str" 
	 */
	public static String removeOuterBrackets(String newick){
		return newick.substring(1,newick.length()-1);
	}
	
	/**
	 * Internal nodes have the following appearance (A:0.5,B:0.5):0.10
	 * @param newick
	 * @return
	 */
	public static boolean isInternalNode(String newick){
		return newick.startsWith("(");
	}
	
	/**
	 * (A:10,(B:20,C:30):50,((E:10,G:10):10,F:10):20 );
	 * -> "A:10", "(B:20,C:30):50" , "((E:10,G:10):10,F:10):20"
	 * @param newick
	 * @return
	 */
	public static List<String> getSubtrees(String newick){
		List<String> subtrees = new ArrayList<String>();
		StringBuilder buffer = new StringBuilder();
		int bracketBalance=0;
		for (int i=0; i<newick.length(); i++){
			char c = newick.charAt(i);
			
			if (c==','){
				if (bracketBalance==0){
					//full subtree
					subtrees.add(buffer.toString());
					buffer.setLength(0); //reset buffer
				} else 
					buffer.append(c);
				 
			} else {
				buffer.append(c);
				if (c=='('){
					bracketBalance++;
				} else if (c==')'){
					bracketBalance--;
				} 
			}
		}
		subtrees.add(buffer.toString()); //rightmost subtree
		
		return subtrees;
	}
	
	/**
	 * (A:10,B:10):20 -> "(A:10,B:10)" , "20"
	 * @param newick
	 * @return
	 */
	public static String [] generateNodeBranchLengthPair(String newick){
		String [] pair = new String[2];
		int indexSemicolon=newick.lastIndexOf(':');
		pair[0]=newick.substring(0,indexSemicolon);
		pair[1]=newick.substring(indexSemicolon+1);
		return pair;
	}
}