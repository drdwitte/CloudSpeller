package be.iminds.cloudspeller.motifalgorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.driver.MotifExtractor;

import be.iminds.cloudspeller.alphabets.CharacterIterator;

import be.iminds.cloudspeller.motifmodels.IUPACMotif;
import be.iminds.cloudspeller.motifmodels.MotifFactory;

import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.phylogenetics.ConservationScoreCalculator;

import be.iminds.cloudspeller.indexing.ISMonkey;
import be.iminds.cloudspeller.indexing.IndexStructure;
import be.iminds.cloudspeller.indexing.IndexStructureFactory;
import be.iminds.cloudspeller.input.Sequence;

/*
Korte Beschrijving discovery algoritme:
Het algoritme heeft de attributen: sequences, index, motifSearchspace, conservationScoreCalculator
en extractor, dewelke voor zich spreken.

in runDiscovery start het eigenlijke algorithme, er wordt een ExactGSMonkey (voor navigatie in
de indexstructuur) aangemaakt die geinitaliseerd staat op de root van de suffixtree, de extractor 
wordt ook gereset. (in het geval van motifContainer is dat een clear() functie)

Vervolgens start een recursieve methode die bestaat uit twee delen: 
exploreSubtree en exploreExtensions dewelke elkaar aanroepen.

In exploreSubtree wordt gecontroleerd het huidige motief zich in de zoekruimte bevindt:
min- en maxLengte, ontaarding, .. dit wordt gedaan in checkSearchSpace en dit levert 4
mogelijkheden op.
Als een motief in de zoekruimte dan wordt een blsScore berekent en wordt getest of die
BLS > BLS_cutoff (survivesBranchAndBound), zoja wordt het opgeslagen in de extractor en worden 
langere motieven onderzocht in exploreExtensions.

In exploreExtensions worden alle motief extensies bekeken dus: mA, mC, mG, mT, mN, enz..
Er wordt telkens een kopie van de monkey genomen (deze bevat alle branchpositions waar een 
ontaard motief matcht) en voor elk van de positie wordt gekeken of het mogelijk is verder te gaan
met een gegeven extension karakter.
Als er matches zijn (hasMatches) dan wordt de nieuwe Monkey naar exploreSubtree gezocht 
waar weerom wordt gecontroleerd of het nieuwe motief binnen de zoekruimte zit, BLS, ...

*getMotifTrail(); -> geef huidige motief

*hasMatches(); -> zijn er matches in de boom?

*grabSuffixes(); -> navigeer van alle branchPositions naar de leaf 
*					en groepeer de suffixes (pattern matching)

*grabInternalNodeInfo(); -> geef de NodeDecorations (=sequenceIDs) van alle knopen 
*							waar het huidige motief matcht

*createClone(); -> deep copy

*jumpTo(Character c); -> ga van huidige positie naar een nieuw (ontaard karakter)

*createDegenerateMonkey(Character degenerateExtension,
			List<ISMonkey> exactMonkeys);
			
			-> ik was aan het proberen de jumpTo te vervangen door iets efficienter maar voorlopig
			is het nog trager. Het idee is om in plaats van met alle 11 de karakters te extenderen
			eerst ACGT te doen en dan gewoon de branchposities van deze 'exactMonkeys'
			samen te voegen, dat bespaart een heleboel jumpTos, voorlopig is het echter nog 
			niet sneller en staat het ergens op de TODO list.


*/

public class DeNovoExactDiscoveryAlgorithm implements DiscoveryAlgorithm {
	
	protected ArrayList<Sequence> sequences;
	protected IndexStructure index;
	protected MotifSearchSpace motifSearchSpace;
	protected ConservationScoreCalculator conservationScoreCalculator;
	protected MotifExtractor extractor;
	//int counter = 0;
	
	public DeNovoExactDiscoveryAlgorithm(ArrayList<Sequence> seqs) {
		this.sequences=seqs;
	}

	@Override
	public void setDataStructure(IndexStructureFactory indexStructureFactory) {
		this.index = indexStructureFactory.createIndexStructure(sequences);
		
	}

	@Override
	public void setSearchSpace(MotifSearchSpace motifSearchSpace) {
		this.motifSearchSpace = motifSearchSpace;
	}

	@Override
	public void setConservationScoreCalculator(ConservationScoreCalculator c) {
		this.conservationScoreCalculator=c;
	}

	@Override
	public void runDiscovery(MotifFactory factory) {
		if (extractor!=null){
			extractor.reset();
			
			exploreSubtreeFast(index.getExactISMonkey(factory,motifSearchSpace.getMaximumDegeneracy()));
		}
	}
	
	public void exploreSubtreeFast(ISMonkey monkey) {

		IUPACMotif currentMotif =(IUPACMotif)monkey.getMotifTrail();
		ConservationScore score = calculateConservationScore(monkey);

		if (!survivesBranchAndBoundCondition(score)){
			return; // regardless of its size, if it does not make the BB condiditions, get out
		}
			
		// only store the motif if it is sufficiently long
		if (currentMotif.length() >= motifSearchSpace.getMinLength()){
			addToMotifContainer(monkey,score);
		}

		if (currentMotif.length() == motifSearchSpace.getMaxLength()){
			return;
		}

		CharacterIterator charIt;
		//check if extensions can be degenerate or not
		if (currentMotif.numberOfDegPositions() == motifSearchSpace.getMaxNumberOfDegeneratePositions()){
			charIt = motifSearchSpace.getAlphabet().exactCharsIterator();
		} else {
			charIt = motifSearchSpace.getAlphabet().getAllCharsIterator();
		}
				
		while (charIt.hasNext()){
			char c = charIt.next();
			monkey.jumpTo(c);
			
			if (monkey.hasMatches())
				exploreSubtreeFast(monkey);
			
			monkey.backtrack();
			
		}

	}
	
	/**
	 * Depth-first de novo discovery algorithm
	 * @param monkey
	 */
	public void exploreSubtree(ISMonkey monkey){
		
		ConservationScore score;
		switch (checkSearchSpace(monkey)){
			
		case NOTINSPACE_NOTEXTENDABLE:
		
			return;

		case NOTINSPACE_EXTENDABLE:
			
			exploreExtensions(monkey);
			return;
			
		case INSPACE_NOTEXTENDABLE:	
			score=calculateConservationScore(monkey);

			if (survivesBranchAndBoundCondition(score)){
				addToMotifContainer(monkey,score);
			}
			return;
			
		case INSPACE_EXTENDABLE:
			score=calculateConservationScore(monkey);

			if (survivesBranchAndBoundCondition(score)){
				addToMotifContainer(monkey,score);
				exploreExtensions(monkey);
			}
			return;
		}
	}
	
	private void addToMotifContainer(ISMonkey monkey,
			ConservationScore score) {
		// Here we take a deep copy because the reference
		// to monkey will be modified as the algorithm progresses
		extractor.add(monkey.getMotifTrail().createDeepCopy(),score);
	}
	
	public int getTempNumberOfMotifs(){
		return extractor.getNumberOfMotifsExtracted();
	}

	/**
	 * Depth first search in tree
	 * @param monkey Iterator alike pattern for index structures 
	 * 
	 */
	
	//@SuppressWarnings({"unused" })
	private void exploreExtensions(ISMonkey monkey) {
	
		
		for (Character c : motifSearchSpace.getAlphabet()){
			ISMonkey newMonkey = monkey.createClone();
			newMonkey.jumpTo(c);
			
			if (newMonkey.hasMatches()){
	
				exploreSubtree(newMonkey);
			}
		}
	}
	
	@SuppressWarnings({"unused" }) //"new Version backup"
	/*private void exploreExtensionsNew(ISMonkey monkey){
		//explore exact extensions
		Map<Character,ISMonkey> exactMonkeyChildren = exploreExactExtensions(monkey);
		
		//explore degenerate extensions
		Alphabet motifAlphabet = motifSearchSpace.getAlphabet();
		List<ISMonkey> exactMatches = new LinkedList<ISMonkey>();
		
		CharacterIterator it = motifAlphabet.degenerateCharsIterator();
		
		while (it.hasNext()){
		
			char degChar=it.next();
			exactMatches.clear(); //avoid allocation
			Iterator<Character> exactCharIt = motifAlphabet.getMatchingCharactersIterator(degChar); 
			while (exactCharIt.hasNext()){
				exactMatches.add(exactMonkeyChildren.get(exactCharIt.next()));
			}
			
			ISMonkey newMonkey = monkey.createDegenerateMonkey(degChar,exactMatches);
			if (newMonkey.hasMatches()){
				exploreSubtree(newMonkey);
			}
		}
		
		//
	}*/

	private Map<Character, ISMonkey> exploreExactExtensions(ISMonkey monkey) {
		Map<Character, ISMonkey> exactExt = new HashMap<Character, ISMonkey>();
		CharacterIterator it = motifSearchSpace.getAlphabet().exactCharsIterator();
		while (it.hasNext()){
			char exChar=it.next();
						
			ISMonkey newMonkey = monkey.createClone();
			newMonkey.jumpTo(exChar);
			
			exactExt.put(exChar,newMonkey);
			
			if (newMonkey.hasMatches()){
				exploreSubtree(newMonkey);
			}
		}
		return exactExt;
	}

	private ConservationScore calculateConservationScore(ISMonkey monkey) {
			return conservationScoreCalculator.calculateScore(monkey.getJointNodeInfo());
	}

	private TrailType checkSearchSpace(ISMonkey monkey) {
		IUPACMotif currentMotif =(IUPACMotif)monkey.getMotifTrail();
		
		if (currentMotif.numberOfDegPositions()>motifSearchSpace.getMaxNumberOfDegeneratePositions()){
			return TrailType.NOTINSPACE_NOTEXTENDABLE;
		} else if (currentMotif.length()<motifSearchSpace.getMinLength()) {
			return TrailType.NOTINSPACE_EXTENDABLE;
		} else if (currentMotif.length()<motifSearchSpace.getMaxLength()){
			return TrailType.INSPACE_EXTENDABLE;
		} else { //k=kmax (can never be bigger)
			return TrailType.INSPACE_NOTEXTENDABLE;
			
		}
	}

	private boolean survivesBranchAndBoundCondition(ConservationScore score){
		return score!=null;
	}

	
	private enum TrailType {
		INSPACE_EXTENDABLE, 
		INSPACE_NOTEXTENDABLE, 
		NOTINSPACE_EXTENDABLE, 
		NOTINSPACE_NOTEXTENDABLE;
	}

	@Override
	public void setMotifExtractor(MotifExtractor extractor) {
		this.extractor=extractor;
	}
}

