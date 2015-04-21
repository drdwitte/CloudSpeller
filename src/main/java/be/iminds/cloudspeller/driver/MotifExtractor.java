package be.iminds.cloudspeller.driver;



import be.iminds.cloudspeller.phylogenetics.ConservationScore;
import be.iminds.cloudspeller.motifmodels.Motif;

public interface MotifExtractor {

	public void add(Motif motif, ConservationScore score);

	public int getNumberOfMotifsExtracted();

	/**
	 * Resets extractor (in case of motif container map is cleared)
	 */
	public void reset();

	public void close();


}

//NOTE idee kan zijn om deels in mapper combining te doen voor de frequentere motieven!
//stel bijvoorbeeld dat we alle 6-meren lokaal aggregeren -> er zal veel overlap zijn tussen
//de verschillende maps, hoe langer de motieven hoe minder kans op overlap dus instant emitter
//lijkt dan het beste, verder kunnen we ook nog een combiner schrijven maar dat zal bij benadering
//een unit functie zijn door de geringe overlap -> partial in-mapper combining lijkt dan eigenlijk
//nog interessanter omdat het het streamen lichtjes kan verminderen.

//het kan ook een idee zijn om ipv group,motief paren te emitten (pairs) om tijdelijk op te slaan en 
// als values motieven te groeperen per permutatiegroep => minder emits (stripes approach)
