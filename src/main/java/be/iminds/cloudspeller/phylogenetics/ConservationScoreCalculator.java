package be.iminds.cloudspeller.phylogenetics;

import be.iminds.cloudspeller.indexing.NodeDecoration;

public interface ConservationScoreCalculator {

	/**
	 * 
	 * @param nodeInfo
	 * @return (conservationScore < cutoffScore)?null:conservationScore;
	 */
	public ConservationScore calculateScore(NodeDecoration nodeInfo);
	public void setCutoff(ConservationScore cutoffScore);
	

}
