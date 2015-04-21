package be.iminds.cloudspeller.indexing;

//SUCCESFUL UNIT TESTS AT 30/5/2013

import be.iminds.cloudspeller.input.BaseSequence;
import be.iminds.cloudspeller.input.Sequence;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import be.iminds.cloudspeller.motifmodels.IUPACFactory;
import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.MMFactory;
import be.iminds.cloudspeller.motifmodels.Motif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;
import org.apache.commons.lang.NotImplementedException;
import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.IUPACAlphabet.IUPACType;

/*
Uitleg bij Generalized SuffixTree klasse:

Deze klasse is vrij omvangrijk: ze bestaat uit een constructiegedeelte en implementeert
een soort iterator-interface (=ISMonkey). De constructie is een stuk gezonder dan in de
C++ code, de ISMonkey laat toe om een erg clean algoritme te schrijven.

Een GST is in een aaneenschakeling van SuffixTreeNodes (lijn 349) die elk een lijst van referenties naar
childnodes bevatten, een STEdge (lijn 329) (incoming edge -> welke karakters?), en afhankelijk
van als het om een leaf gaat een lijst van Suffixes of als het om een interne node gaat een 
NodeDecoration (sequentiesIndices nu geimplementeerd als een bitset)

De constructie wordt volledig afgehandeld in de private class ConstructionAlgorithm (lijn 133),
het algoritme is anders (=beter) dan in de originele code. Het probleem met de originele code
was dat je na de sorteeroperatie nog niet weet welke knoop overeenstemt met een bepaalde SuffixBucket
Nu wordt er steeds gesorteerd tot op een split, vervolgens wordt voor elke SuffixBucket gekeken
hoe lang de tak in de boom moet verlengd worden tot de volgende split.
Het resultaat van dit algoritme is getest met de toString() methode die de boom wegschrijft in 
dotFormat hetgeen met commandline: dot -Tpdf testGraph.txt > testGraph.pdf kan worden gevisualiseerd.

Implementatie van ISMonkey (lijn 447)
Deze navigatie-klasse dient in de eerste plaats als container voor de matches van een gegeven 
motief (-> BranchPosition (lijn 512) en een Motif Trail (dit is het pad door de zoekruimte = huidig
ontaard motief).
De ISMonkey vereist implementatie van volgende functies:

*/
public class GeneralizedSuffixTree implements IndexStructure {
	
	public static final String SENTINEL = "s";
	
	private SuffixTreeNode root;
	private NodeDecorationFactory nodeDecorationFactory;
	private ArrayList<Sequence> text;
	private Map<Character,Integer> charChildMap = new HashMap<Character,Integer>();

	private int numberOfForwardSequences;	
	
	/**
	 * 
	 * @return Sequences with sentinel character + reverse complements added+
	 * @return Replace non alphabet chars with sentinels
	 */
	public static ArrayList<Sequence> preprocessSequences(ArrayList<Sequence> sequences, boolean withReverseComplements) {
	
		ArrayList<Sequence> prepSeqs = new ArrayList<Sequence>();
		
		int numberOfForwardSequences=sequences.size();
		for (int i = 0; i < numberOfForwardSequences; i++) {
			prepSeqs.add(sequences.get(i).getClone());
		}
		
						
		if (withReverseComplements){ //reverse complement implies baseSeq
			
			for (int i = 0; i < numberOfForwardSequences; i++) {
				Sequence currentSeq = sequences.get(i);
				BaseSequence baseSeq = (BaseSequence) currentSeq;
				prepSeqs.add(baseSeq.generateReverseComplement());
			}
		}
		
	
		addSentinels(prepSeqs);

		return prepSeqs;
	}
	
	private static void addSentinels(ArrayList<Sequence> sequences) {

		String validChars = sequences.get(0).getSequenceAlphabet().getAllChars()+SENTINEL;
		
		for (int i = 0; i < sequences.size(); i++) {

			Sequence currentSeq = sequences.get(i);
			//replace non alphabet characters
			for (int j=0; j<currentSeq.length(); j++){
				char c = currentSeq.charAt(j);
				
				if (Character.isLowerCase(c)){
					c = Character.toUpperCase(c);
					currentSeq.replace(j,j+1,""+c);
				}
				
				if (validChars.indexOf(c)==-1){
					currentSeq.replace(j,j+1,SENTINEL);
				} else {}
			}
			
			currentSeq.append(SENTINEL);
		}
	}

	/**
	 * Generalized Suffix Tree Constructor: run WOTD algorithm and add NodeDecoration 
	 */
	public GeneralizedSuffixTree(ArrayList<Sequence> sequences, boolean withReverseComplements, int maxDepth, NodeDecorationFactory decorationFac) {
		
		this.nodeDecorationFactory=decorationFac;
		
		int i = 0;
		for (Character c : sequences.get(0).getSequenceAlphabet()) {
			charChildMap.put(c, i++);
		}
		charChildMap.put(SENTINEL.charAt(0),i);
		this.text = preprocessSequences(sequences, withReverseComplements);
		this.numberOfForwardSequences= withReverseComplements?text.size()/2:text.size();
		
		// run WOTD algorithm
		List<Suffix> suffixes = generateSuffixes();
		root = new InternalNode(new STEdge(""),suffixes);
		ConstructionAlgorithm constructionAlgorithm = new ConstructionAlgorithm(
				root, maxDepth);
		constructionAlgorithm.runWOTDAlgorithm(suffixes);

	}
	
	/**
	 * Generalized Suffix Tree Constructor: run WOTD algorithm and add NodeDecoration
	 * for a user defined set of suffixes (can be used for for be.iminds.cloudspeller.indexing in aligned sequences
	 * NOTE: currently no reverse strands supported
	 */
	public GeneralizedSuffixTree(ArrayList<Sequence> preprocessedSequences, List<Suffix> suffixes, int maxDepth, NodeDecorationFactory decorationFac) {
	
		this.nodeDecorationFactory=decorationFac;
		
		int i = 0;
		for (Character c : preprocessedSequences.get(0).getSequenceAlphabet()) {
			charChildMap.put(c, i++);
		}
		charChildMap.put(SENTINEL.charAt(0),i);
		this.text = preprocessedSequences;
		this.numberOfForwardSequences= text.size();
		
		// run WOTD algorithm
		root = new InternalNode(new STEdge(""),suffixes);
		ConstructionAlgorithm constructionAlgorithm = new ConstructionAlgorithm(
				root, maxDepth);
		constructionAlgorithm.runWOTDAlgorithm(suffixes);
		
	}
	
	
	
	public List<Suffix> generateSuffixes(){
		List<Suffix> suffixes = new ArrayList<Suffix>(); 
		for (int i=0; i<text.size(); i++){
			for (int j=0; j<text.get(i).length(); j++){
				suffixes.add(new Suffix(i,j));
			}
		}
		return suffixes;
	}
		
	private class ConstructionAlgorithm{
		private SuffixTreeNode root;
		private int maxDepth;
		
		
		public ConstructionAlgorithm(SuffixTreeNode root, int maxDepth) {
			this.root = root;
			this.maxDepth = maxDepth;
		}
		
		/**
		 * Write-Only Top Down algorithm adapted from Kurz (original constr. algorithm is lazy and for ST)
		 */
		public void runWOTDAlgorithm(List<Suffix> suffixes) {
			SuffixBucket initialBucket = new SuffixBucket(suffixes,"");
			createSubtree(root,0,initialBucket.sortBucket(0)); //note: initial minimal 2 buckets: sentinelchar and other char
			
			//TODO root = leaf case is niet supported!! (lege data enkel sentinels!)
			
		}
		
		
		private class SuffixBucket {
			

			private String bucketPrefix = "";
			private List<Suffix> suffixList;
			
			public SuffixBucket(List<Suffix> suffixList, String bucketPrefix){
				this.suffixList=suffixList;
				setPrefix(bucketPrefix);
			}

			public List<SuffixBucket> sortBucket(int depth) {
				SuffixBucket [] bucketArray = new SuffixBucket[charChildMap.size()];
				
				for (Entry<Character, Integer> entry : charChildMap.entrySet()){
					bucketArray[entry.getValue()] = new SuffixBucket(new ArrayList<Suffix>(),""+entry.getKey());
				}
							
				for (Suffix suffix : suffixList ){
					Sequence sequence = text.get(suffix.getSequenceID());
					int position = suffix.getSequencePosition() + depth;
					if (position>=sequence.length()){
						continue;
					}
					char currentChar = sequence.charAt(position);
					
					bucketArray[charChildMap.get(currentChar)].addSuffix(suffix);
				}
				
				//remove empty buckets
				List<SuffixBucket> buckets = new ArrayList<SuffixBucket>();
				for (SuffixBucket bucket : bucketArray){
					if (bucket.getSuffixBucket().size()!=0){
						buckets.add(bucket);
					}
				}
				
				
				
				if (buckets.size() == 1){ //nothing sorted extend bucket prefix
					String extensionChar = buckets.get(0).getPrefix();
					buckets.get(0).setPrefix(bucketPrefix+extensionChar);
				}
				
				return buckets;
			}

			private void setPrefix(String string) {
				this.bucketPrefix=string;
				
			}

			private void addSuffix(Suffix suffix) {
				suffixList.add(suffix);
				
			}

			public String getPrefix() {
				return bucketPrefix;
			}
			

			
			public List<Suffix> getSuffixBucket(){
				return suffixList;
			}

			public boolean isSentinelBucket() {
				return bucketPrefix.endsWith(SENTINEL);
			}
			
			@Override
			public String toString() {
				return bucketPrefix+"\n"+suffixList;
			}
	
		}
		
		private Leaf createLeafFromBucket(SuffixBucket bucket){
			STEdge edge = new STEdge(bucket.getPrefix());
			Leaf leaf = new Leaf(edge,bucket.getSuffixBucket());
			return leaf;
		}
		
		private InternalNode createInternalNodeFromBucket(SuffixBucket bucket) {
			STEdge edge = new STEdge(bucket.getPrefix());
			InternalNode iNode = new InternalNode(edge,bucket.getSuffixBucket());
			return iNode;
		}
		
		private void connectParentAndChild(SuffixTreeNode parent, SuffixTreeNode child) {
			InternalNode parentInternal = (InternalNode) parent;
			parentInternal.connectChild(child.charAtEdge(0), child);
		}
		
		
		/**
		 * 
		 * @param parent
		 * @param parentStringDepth current position of suffix is Sj+stringDepth-1
		 * => initial children have stringdepth 1 => sort position 0
		 * @param buckets are sorted at stringdepth 1 => position 0
		 */
		private void createSubtree(SuffixTreeNode parent, int parentStringDepth, List<SuffixBucket> buckets){
			
			int currentChildDepth=parentStringDepth+1;
			
			if (currentChildDepth == maxDepth) { // create all leafs

				for (SuffixBucket bucket : buckets) {
					Leaf leaf = createLeafFromBucket(bucket);
					connectParentAndChild(parent,leaf);
				}
				
			} else {

				OuterForLoop:for (SuffixBucket bucket : buckets) {
					currentChildDepth=parentStringDepth+1;
					int charIDToSort=parentStringDepth+1; //buckets correspond to the sorting of 
														   //suffixes one char below parent	
					
					if (bucket.isSentinelBucket()) {

						Leaf sentinelLeaf = createLeafFromBucket(bucket);

						connectParentAndChild(parent,sentinelLeaf);
						continue OuterForLoop;
					}

					List<SuffixBucket> newBuckets;

					newBuckets = bucket.sortBucket(charIDToSort++);

					SuffixBucket singleBucket = bucket;
					while (newBuckets.size() == 1) {

						currentChildDepth++;
						singleBucket = newBuckets.get(0);
						
						if (currentChildDepth == maxDepth || singleBucket.isSentinelBucket()) {

							Leaf leaf = createLeafFromBucket(singleBucket);
							connectParentAndChild(parent, leaf);
							continue OuterForLoop; // goto next bucket
						}

						newBuckets = singleBucket.sortBucket(charIDToSort++);


					}
					
					InternalNode internalNode = createInternalNodeFromBucket(singleBucket);
					connectParentAndChild(parent, internalNode);
					createSubtree(internalNode, currentChildDepth , newBuckets);
				}
			}
		}
	
	}

	//Substitute if beneficial: currently strings on edges
	/*private abstract class SuffixTreeNode {
		protected short startID;
		protected short stopID;
		protected byte seqID;
		
		public SuffixTreeNode(short startID, short stopID, byte seqID) {
			this.startID = startID;
			this.stopID = stopID;
			this.seqID = seqID;
		}

		public abstract boolean isLeaf();

		public Character charAt(int i) {
			return text.get(seqID).charAt(startID+i);
			
		}
	}*/
	
	private class STEdge {
		private String label;
		
		public STEdge(String label){
			this.label=label;
		}
		
		public char charAt(int i){
			return label.charAt(i);
		}
		
		public int getLength(){
			return label.length();
		}
		
		public String getEdgeLabel(){
			return label;
		}
	}
	
	private abstract class SuffixTreeNode {
		protected STEdge incomingEdge;
		protected NodeDecoration nodeInfo;
		
		public SuffixTreeNode(STEdge incomingEdge, List<Suffix> suffixes) {
			this.incomingEdge = incomingEdge;
			generateNodeDecoration(suffixes);
		}

		public abstract boolean isLeaf();

		public Character charAtEdge(int i) {
			return incomingEdge.charAt(i);
		}
		
		public STEdge getEdge(){
			return incomingEdge;
		}
		
		private void generateNodeDecoration(List<Suffix> suffixes){
			nodeInfo = nodeDecorationFactory.createNodeDecoration(numberOfForwardSequences);
			nodeInfo.processSuffixes(suffixes);
		}

		public NodeDecoration getNodeInfo() {
			return nodeInfo;
		}

		@Override
		public String toString() {
			return "("+incomingEdge.label+")["+nodeInfo+"] "+this.hashCode();
		}
		
		
	}
	
	private class InternalNode extends SuffixTreeNode {
		private SuffixTreeNode [] children;
		
		public InternalNode(STEdge edge, List<Suffix> suffixes){
			super(edge,suffixes);
			children = new SuffixTreeNode[charChildMap.size()];                     
		}
		
		public void connectChild(Character c, SuffixTreeNode child){
			children[charChildMap.get(c)]=child;
		}

		@Override
		public boolean isLeaf() {
			return false;
		}
		
		public SuffixTreeNode getChild(Character c){
			return children[charChildMap.get(c)];
		}
		
		
	}
	
	private class Leaf extends SuffixTreeNode {
		
		private List<Suffix> suffixes = new ArrayList<Suffix>();
		
		public Leaf(STEdge edge, List<Suffix> suffixes) {
			super(edge,suffixes);
			addSuffixes(suffixes);
		}

		private void addSuffixes(List<Suffix> suffixList) {
			suffixes.addAll(suffixList);
		}

		@Override
		public boolean isLeaf() {
			return true;
		}
		
		public List<Suffix> getSuffixes(){
			return suffixes;
		}
		
	}
	
	protected class BranchPosition {
		
		private static final int BRANCHPOSTION_ROOT=-1;
		
		private SuffixTreeNode nodeBelow;
		private int edgePosition;
		
		public BranchPosition() {
			nodeBelow = null;
			edgePosition = 0;
		}
		
		public BranchPosition(SuffixTreeNode nodeBelow, int edgePosition){
			setNodeBelow(nodeBelow);
			setEdgePosition(edgePosition);
		}
		
		
		public STEdge getEdge() {
			return nodeBelow.incomingEdge;
		}


		public SuffixTreeNode getNodeBelow() {
			return nodeBelow;
		}

		/**
		 * @return true if postion points at node (which means last character on edge
		 */
		public boolean inNode() {
			return (edgePosition+1) == nodeBelow.incomingEdge.getLength();
		}

		public void setNodeBelow(SuffixTreeNode nodeBelow) {
			this.nodeBelow = nodeBelow;
		}

		public int getEdgePosition() {
			return edgePosition;
		}

		public void setEdgePosition(int edgePosition) {
			this.edgePosition = edgePosition;
		}
		
		public void setData(SuffixTreeNode nodeBelow, int edgePosition) {
			this.nodeBelow = nodeBelow;
			this.edgePosition = edgePosition;
		}
		
		public String toString(){
			return nodeBelow.toString()+"\t"+edgePosition;
		}
		
	}
	
	private abstract class GSMonkey implements ISMonkey {
		
		protected List<BranchPosition> positions = new LinkedList<BranchPosition>();
		protected Motif trail;
		protected Alphabet motifAlphabet;
		
		public GSMonkey(Motif trail){
			this.trail = trail;
		}
		
		@Override
		public ArrayList<NodeDecoration> grabInternalNodeInfo() {
			ArrayList<NodeDecoration> nodeInfo = new ArrayList<NodeDecoration>();
			for (BranchPosition pos : positions){
				nodeInfo.add(pos.getNodeBelow().getNodeInfo());
			}
			return nodeInfo;
		}
		
		@Override
		public NodeDecoration getJointNodeInfo() {
			NodeDecoration jointNodeInfo = positions.get(0).getNodeBelow().getNodeInfo().createClone();
			
			for (BranchPosition pos : positions)
				jointNodeInfo.joinWith(pos.getNodeBelow().getNodeInfo());				
			
			return jointNodeInfo;
		}

		@Override
		public List<Suffix> grabSuffixes() { 
			List<Suffix> allSuffixes = new ArrayList<Suffix>();
			for (BranchPosition position : positions ){
				
				grabSuffixes(position.getNodeBelow(),allSuffixes);
				
			}
			return allSuffixes;
		}
		
		private void grabSuffixes(SuffixTreeNode node, List<Suffix> allSuffixes){
			
			if (node.isLeaf()){
				Leaf leaf = (Leaf) node;
				allSuffixes.addAll(leaf.getSuffixes());
				return;
			} else {
				InternalNode iNode = (InternalNode) node;
				for (SuffixTreeNode child : iNode.children){
					
					if (child!=null){
						
						grabSuffixes(child,allSuffixes);
					}
				}
			}
			return;
		}
		
		@Override
		public boolean hasMatches() {
			return positions.size()>0;
		}

		@Override
		public abstract void jumpTo(Character c);

		@Override
		public Motif getMotifTrail() {
			return trail;
		}
		
		@Override
		public abstract ISMonkey createClone(); 
		
	}
	
	private class ExactGSMonkey extends GSMonkey {

		/**
		 * Constructor for ExactGSMonkey initializing it to the root of the GST
		 * @param factory Motif factory used to initialized the motif trail and the motif alphabet
		 */
		public ExactGSMonkey(MotifFactory factory) {
			super(factory.createEmptyMotif());
			positions.add(new BranchPosition(root,BranchPosition.BRANCHPOSTION_ROOT));
			this.motifAlphabet=factory.getAlphabet();
		}
		
		private ExactGSMonkey(){
			super(null);
		}

		public void jumpTo(Character c) {
						
			List<BranchPosition> newPositions = new LinkedList<BranchPosition>();
			Iterator<Character> extensionIterator=motifAlphabet.getMatchingCharactersIterator(c);
			
			while (extensionIterator.hasNext()) {
				Character extension = extensionIterator.next();

				// update positions
				for (BranchPosition pos : positions) {

					// check if position is in node or inside edge
					if (pos.inNode()) {
						

						if (pos.getNodeBelow().isLeaf()) { // cannot extend past leaf
							continue;
						} else {
							InternalNode parentNode = (InternalNode) pos.getNodeBelow();
							SuffixTreeNode child = parentNode.getChild(extension);
							
							if (child != null) { // check if branch exists to extension character
								newPositions.add(new BranchPosition(child, 0));
								
							} else {}
						}

					} else {
						
						Character nextCharOnEdge = pos.getEdge().charAt(pos.edgePosition + 1);
						if (nextCharOnEdge.equals(extension)) {
							newPositions.add(new BranchPosition(pos.getNodeBelow(),pos.getEdgePosition() + 1));
						} else {}
					}
				}
			}

			// update trail
			trail.append(c);
			this.positions=newPositions;
		}

		@Override
		public ISMonkey createClone() {
			ExactGSMonkey gsMonkey = new ExactGSMonkey();
			gsMonkey.positions =  new LinkedList<BranchPosition>(this.positions);
			gsMonkey.trail = this.trail.createDeepCopy();
			gsMonkey.motifAlphabet= this.motifAlphabet;
			return gsMonkey;
		}

		/*@Override
		public ISMonkey createDegenerateMonkey(Character degenerateExtension,
				List<ISMonkey> exactMonkeys) {
			ExactGSMonkey newMonkey = new ExactGSMonkey();
			newMonkey.positions = new LinkedList<BranchPosition>();
			for (ISMonkey m : exactMonkeys){
				newMonkey.positions.addAll(((GSMonkey)m).positions);
			}
			
			newMonkey.trail = this.trail.createDeepCopy();
			newMonkey.trail.append(degenerateExtension);
			
			newMonkey.motifAlphabet=this.motifAlphabet;
			
			return newMonkey;
		}*/

		@Override
		public void backtrack() {
			throw new UnsupportedOperationException();
		}	
	}
	
	private class GSMonkeyFast implements ISMonkey {

		protected BranchPosition[][] positions;
		protected int[] posSize;

		protected Motif trail;
		protected Alphabet motifAlphabet;
		protected int currentLevel;

		public GSMonkeyFast(MotifFactory factory, int maxNavigationDepth, int maxNumberOfMatchingSites) {
		
			trail = factory.createEmptyMotif();
			positions = new BranchPosition[maxNavigationDepth][]; 
			posSize = new int[maxNavigationDepth]; 

			for (int i = 0; i < maxNavigationDepth; i++) {
				positions[i] = new BranchPosition[maxNumberOfMatchingSites]; 
				for (int j = 0; j < maxNumberOfMatchingSites; j++)
					positions[i][j] = new BranchPosition();				
				posSize[i] = 0;
			}

			positions[0][0].setData(root, BranchPosition.BRANCHPOSTION_ROOT);
			posSize[0] = 1;

			currentLevel = 0;

			motifAlphabet=factory.getAlphabet();
		}
	
		@Override
		public Motif getMotifTrail() {
			return trail;
		}

		@Override
		public boolean hasMatches() {
			return posSize[currentLevel] > 0;
		}

		@Override
		public void jumpTo(Character c) {
			
			Iterator<Character> extensionIterator=motifAlphabet.getMatchingCharactersIterator(c);
			while (extensionIterator.hasNext()) {
				Character extension = extensionIterator.next();

				// update positions
				for (int i = 0; i < posSize[currentLevel]; i++) {
					BranchPosition pos = positions[currentLevel][i];

					// check if position is in node or inside edge
					if (pos.inNode()) {
						if (pos.getNodeBelow().isLeaf()) { // cannot extend past leaf
							continue;
						} else {
							InternalNode parentNode = (InternalNode) pos.getNodeBelow();
							SuffixTreeNode child = parentNode.getChild(extension);
							
							if (child != null) { // check if branch exists to extension character
								positions[currentLevel + 1][posSize[currentLevel + 1]].setData(child, 0);
								posSize[currentLevel + 1]++;
							}
						}
					} else {						
						Character nextCharOnEdge = pos.getEdge().charAt(pos.edgePosition + 1);
						if (nextCharOnEdge.equals(extension)) {
							positions[currentLevel + 1][posSize[currentLevel + 1]].setData(pos.getNodeBelow(),pos.getEdgePosition() + 1);
							posSize[currentLevel + 1]++;
						}
					}
				}
			}
			// update trail
			trail.append(c);
			currentLevel++;
		}

		/*@Override
		public ISMonkey createDegenerateMonkey(Character degenerateExtension,
				List<ISMonkey> exactMonkeys) {
			throw new NotImplementedException();
		}*/

		@Override
		public void backtrack() {
			posSize[currentLevel] = 0;
			currentLevel--;
			trail.pop();
		}
		
		@Override
		public ArrayList<NodeDecoration> grabInternalNodeInfo() {
			ArrayList<NodeDecoration> nodeInfo = new ArrayList<NodeDecoration>();
			for (int i = 0; i < posSize[currentLevel]; i++)
				nodeInfo.add(positions[currentLevel][i].getNodeBelow().getNodeInfo());
			return nodeInfo;
		}

		@Override
		public List<Suffix> grabSuffixes() { 
			List<Suffix> allSuffixes = new ArrayList<Suffix>();
			for (int i = 0; i < posSize[currentLevel]; i++)				
				grabSuffixes(positions[currentLevel][i].getNodeBelow(),allSuffixes);				
			return allSuffixes;
		}
		
		@Override
		public NodeDecoration getJointNodeInfo() {					
			NodeDecoration jointNodeInfo = positions[currentLevel][0].getNodeBelow().getNodeInfo().createClone();
			
			for (int i = 1; i< posSize[currentLevel]; i++)
				jointNodeInfo.joinWith(positions[currentLevel][i].getNodeBelow().getNodeInfo());
			
			return jointNodeInfo;
		}
		
		private void grabSuffixes(SuffixTreeNode node, List<Suffix> allSuffixes){
			
			if (node.isLeaf()){
				Leaf leaf = (Leaf) node;
				allSuffixes.addAll(leaf.getSuffixes());
				return;
			} else {
				InternalNode iNode = (InternalNode) node;
				for (SuffixTreeNode child : iNode.children){
					
					if (child!=null){
						
						grabSuffixes(child,allSuffixes);
					}
				}
			}
			return;
		}

		@Override
		public ISMonkey createClone() {
			throw new UnsupportedOperationException();
		}
	}
	
	private class InexactGSMonkey extends GSMonkey {

		public InexactGSMonkey() {
			super(new MMFactory().createEmptyMotif());
			throw new NotImplementedException();
		}

		@Override
		public void jumpTo(Character c) {
			throw new NotImplementedException();
		}

		@Override
		public ISMonkey createClone() {
			throw new NotImplementedException();
		}

		/*@Override
		public ISMonkey createDegenerateMonkey(Character degenerateExtension,
				List<ISMonkey> exactMonkeys) {
			throw new UnsupportedOperationException();
		}*/

		@Override
		public void backtrack() {
			throw new UnsupportedOperationException();			
		}
	}

	@Override
	public ISMonkey getExactISMonkey(MotifFactory motifFactory, int maxDegeneracy) {
	    return new GSMonkeyFast(motifFactory,motifFactory.getMaxLength()+1,maxDegeneracy);
		//return new ExactGSMonkey(motifFactory);
	}

	@Override
	public ISMonkey getInexactISMonkey() {
		return new InexactGSMonkey();
	}

	@Override
	/**
	 * Generate dot representation of GST
	 * The visual reprenstation can be generated on the commandline: dot -Tpdf testGraph.txt > testGraph.pdf
	 * (part of graphviz package) 
	 */
	public String toString() {
		StringBuilder s = new StringBuilder();
		String start = "digraph G {";
		String nl = "\n";
		String stop = "}";
			
			
		s.append(start); s.append(nl);
		
		List<SuffixTreeNode> activeNodes = new LinkedList<SuffixTreeNode>();
		activeNodes.add(root);
		while (activeNodes.size()>0){
			SuffixTreeNode currentNode = activeNodes.get(0);
			activeNodes.remove(0);
			
			if (currentNode.isLeaf()){
				continue;
			} else { //connect to children
				InternalNode iNode = (InternalNode) currentNode;
				
				for (Map.Entry<Character, Integer> entry : charChildMap.entrySet()){
					SuffixTreeNode child = iNode.getChild(entry.getKey()); 
					if (child!=null){
						//draw edge
						s.append(drawEdge(iNode,child)); s.append(nl); 
						activeNodes.add(child); 
					}
				}
			}
		}		
		
		s.append(stop); s.append(nl);		
		return s.toString();		
	}

	/**
	 * Generate dot representation for a suffix tree edge: 
	 * A1 -> A2 [label=f];
	 */
	private String drawEdge(InternalNode iNode, SuffixTreeNode child) {
		String e = " -> ";
		StringBuilder s= new StringBuilder();
		//label
		s.append(drawNode(iNode));
		s.append(e);
		s.append(drawNode(child));
		s.append(" [label=");
		s.append(child.getEdge().getEdgeLabel());
		s.append("];");
		
		return s.toString();
	}
	
	private String drawNode(SuffixTreeNode n){
		StringBuilder s= new StringBuilder();
		if (n.isLeaf()){
			s.append("l");
		} else {
			s.append("n");
		}
		s.append(n.hashCode());
		s.append("_");
		s.append(n.getNodeInfo().toString());
		
		return s.toString();
	}

	
	@Override
	public List<Suffix> matchExactPattern(IUPACMotif pattern) {
		
		ISMonkey monkey = new ExactGSMonkey(new IUPACFactory(IUPACType.FULL));

		for (int i=0; i<pattern.length(); i++){
			monkey.jumpTo(pattern.charAt(i));
			
			if (!monkey.hasMatches()){
				return null;
			}
		}
		
		return monkey.grabSuffixes();

	}


	public ArrayList<Sequence> getPreprocessedSequences() {
		return text;
	}
	

	
}


	
