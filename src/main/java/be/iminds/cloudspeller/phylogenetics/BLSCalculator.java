package be.iminds.cloudspeller.phylogenetics;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.lang.NotImplementedException;

import be.iminds.cloudspeller.indexing.BitSetDecoration;
import be.iminds.cloudspeller.indexing.NodeDecoration;
import be.iminds.cloudspeller.input.Gene;
import be.iminds.cloudspeller.input.GeneFamily;

public class BLSCalculator implements ConservationScoreCalculator {

	private ConservationScore cutoffScore = new BLS(0);
	BLS precalculatedMap [] = null;
	
	/**
	 * Constructor
	 * @param family
	 */
	public BLSCalculator(GeneFamily family){
		
		int numberOfBitsInInt=32;
		int numberOfGenes = family.getGenes().size();
		if (numberOfGenes>numberOfBitsInInt){
			System.out.println(family.getFamilyName());
			System.out.println(family.getGenes());
			throw new NotImplementedException("Only support up until 32 genes per orthologous genefamily");
		}
				
		String substSeqIDs=generateModifiedNewick(family);
		SeqIDTree seqIDTree = new SeqIDTree(substSeqIDs,family.getNumberOfGenes());
		seqIDTree.renormalize();
		precalculateScores(numberOfGenes, seqIDTree);
		
	}
	
	//SETTERS
	@Override
	public void setCutoff(ConservationScore cutoffScore) {
		this.cutoffScore = cutoffScore;
		
	}
	
	//GETTERS
	
	
	//METHODS
	/**
	 * Substitute genes with sequence identifiers (integers)
	 */
	public static String generateModifiedNewick(GeneFamily gf) {
		ArrayList<Gene> genes=gf.getGenes();
		String modifiedString=gf.getNewick().toString();
		for (int i=0; i<genes.size(); i++){
			modifiedString=modifiedString.replaceFirst(genes.get(i).getID(), String.valueOf(i));
		}
	
		return modifiedString;
	}
	
	@Override
	public ConservationScore calculateScore(NodeDecoration nodeInfo) {
		
		ConservationScore score = precalculatedMap[nodeInfo.toIntegerBitRepresentation()];
		
		if (score.compareTo(cutoffScore)<0){
			return null;
		} else 
			return score;
	}
	
	private void precalculateScores(int numberOfGenes, SeqIDTree seqIDTree){
		int nCombinations=(int)Math.pow(2,numberOfGenes);
		precalculatedMap = new BLS[nCombinations];
		precalculatedMap[0]=new BLS(0);
		for (int i=1; i<nCombinations; i++){
			Set<Integer> IDs = new BitSetDecoration(numberOfGenes,i).getIDs();
			precalculatedMap[i] = (BLS) seqIDTree.calculateScore(IDs);
		}
	}

	/*
	private ConservationScore generateRandomScore() {
		BLS score=new BLS((int) (Math.random()*101));
		
		if (score.compareTo(cutoffScore)<0){
			return null;
		} else {
			return score;
		}
	}*/
	
	@Override
	public String toString() {
		
		return "BLSCalculator"; // seqIDTree.toString();
	}
	
	private class SeqIDTree {
		
		PhyloNode [] leafs;
		
		/**
		 * @param seqIDNewickString newick string with gene 
		 * names substituted with sequenceIDs
		 */
		public SeqIDTree(String newickString, int numSeq){
			leafs = new PhyloNode[numSeq];
			int lastID=newickString.length()-1; //remove final ;
			String newickRoot=NewickParser.removeOuterBrackets(newickString.substring(0,lastID));
			expandInternalNode(createRoot(newickRoot));
		}
		
		public void renormalize() {
			ConservationScore normalizationFactor = calculateSumOfEdges();
			
			double doubleValue =  ((BLS)normalizationFactor).doubleValue()/100;
			renormalizeNodes(doubleValue);
			
		}
		
		private void renormalizeNodes(double normalizationFactor){
			
			//initialize active nodes
			SortedSet<PhyloNode> activeNodes = new TreeSet<PhyloNode>();
			for (PhyloNode leaf : leafs){
				activeNodes.add(leaf);
			}
			
			while (activeNodes.size()>1){
				PhyloNode deepestNode=activeNodes.first();
				double newEdgeLength = deepestNode.getEdgeLength() / normalizationFactor;
				deepestNode.setEdgeLength(newEdgeLength);				
				activeNodes.add(deepestNode.getParent());
				activeNodes.remove(deepestNode);
			}
			
		}
		

		private void expandInternalNode(PhyloNode parent) {
			List<String> subtrees=NewickParser.getSubtrees(parent.getNewick());
			for (String subtree : subtrees){
				if (NewickParser.isInternalNode(subtree)){ 
					//create internal node
					String [] pair=NewickParser.generateNodeBranchLengthPair(subtree);
					double edgeLength=Double.parseDouble(pair[1]);
					String nodeID=NewickParser.removeOuterBrackets(pair[0]);
					expandInternalNode(new PhyloNode(parent, edgeLength,nodeID));
				} else { 
					//create leaf
					String [] pair=NewickParser.generateNodeBranchLengthPair(subtree);
					int seqID=Integer.parseInt(pair[0]);
					double edgeLength=Double.parseDouble(pair[1]);
					leafs[seqID]=new PhyloNode(parent,edgeLength,pair[0]);
				}
			}
			
			
		} //END SeqIDTree::expandInternalNode

		/**
		 * Warning changed BLSsum definition -> common edge is always added but can be easily
		 * modified! (for example remaining common edge is counted *1/2?
		 * @param sequenceIDs
		 * @return
		 * NOTE: sequenceIDSet must be >0 otherwise crashes (should never occur!)
		 */
		public ConservationScore calculateScore(Set<Integer> sequenceIDs){
			
			//initialize active nodes
			SortedSet<PhyloNode> activeNodes = new TreeSet<PhyloNode>();

			for (Integer i : sequenceIDs){
				activeNodes.add(leafs[i]);
			}
			
			double sumOfBranchLengths= processActiveNodes(activeNodes);
						
			return new BLS((int) Math.round(100*sumOfBranchLengths));
		}
		
		public ConservationScore calculateSumOfEdges(){
			//initialize active nodes
			SortedSet<PhyloNode> activeNodes = new TreeSet<PhyloNode>();
			for (PhyloNode leaf : leafs){
				activeNodes.add(leaf);
			}
			return new BLS((int) Math.round(100*processActiveNodes(activeNodes)));
		}
		
		private double processActiveNodes(SortedSet<PhyloNode> activeNodes){
			double sumOfBranchLengths=0;
						
			while (activeNodes.size()>1){
				PhyloNode deepestNode=activeNodes.first();
				sumOfBranchLengths+=deepestNode.getEdgeLength();
				activeNodes.add(deepestNode.getParent());
				activeNodes.remove(deepestNode);
			}
			//sumOfBranchLengths+=activeNodes.first().getEdgeLength(); //no ancestor count
			return sumOfBranchLengths;
		}
		
		public PhyloNode createRoot(String newickSubtree){
			return new PhyloNode(newickSubtree);
		}
		
		@Override 
		public String toString(){
			StringBuilder sB = new StringBuilder();
			String arrow="->";
			for (int i=0; i<leafs.length; i++){
				leafs[i].getEdgeLength();
				sB.append("leaf"+i);
				sB.append(arrow);
				sB.append("(");
				sB.append(leafs[i].getEdgeLength());
				sB.append(")");
				PhyloNode parent = leafs[i].getParent();
				
				while (parent.getParent()!=null){
					sB.append(arrow);
					sB.append("(");
					sB.append(parent.getEdgeLength());
					sB.append(")");
					parent=parent.getParent();
				}
				sB.append("R");
				sB.append("\n");
				
				
			}
			return sB.toString();
		}
		
		
		class PhyloNode implements Comparable<PhyloNode> {
			
			private PhyloNode parent;
			private double nodeDepth;
			private double edgeLength;
			private String newickSubtree;
			
			/**
			 * Create Root node
			 */
			public PhyloNode(String newickSubtree){
				this.newickSubtree=newickSubtree;
			}
			
			public void setEdgeLength(double newEdgeLength) {
				this.edgeLength = newEdgeLength;
			}

			/**
			 * Create internal/leaf node in SeqIDTree
			 * @param parent Link to parent node
			 * @param edgeLength Length of incoming edge
			 */
			public PhyloNode(PhyloNode parent, double edgeLength, String newickSubtree){
				this.parent=parent;
				this.nodeDepth=parent.nodeDepth+edgeLength;
				this.edgeLength=edgeLength;
				this.newickSubtree=newickSubtree;
			}
			
			//GETTERS
			public String getNewick() {
				return newickSubtree;
			}
		
			
			public PhyloNode getParent() {
				return parent;
			}

			public double getEdgeLength() {
				return edgeLength;
			}

			//METHODS
			@Override 
			/**
			 * Node with biggest depth is below (=behind) other node -> negative
			 */
			public int compareTo(PhyloNode o) {
				if (nodeDepth<o.nodeDepth)
					return +1; 
				else if (nodeDepth>o.nodeDepth){
					return -1;
				} else 
					return newickSubtree.compareTo(o.newickSubtree);
					
			}
			
		}//END SeqIDTree::Phylonode
		
	}//END SeqIDTree


}//END BLSCalculator


