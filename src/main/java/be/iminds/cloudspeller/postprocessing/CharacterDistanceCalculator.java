package be.iminds.cloudspeller.postprocessing;

public interface CharacterDistanceCalculator {

	public double calculateDistance(char c1, char c2);
	public char getBulkIndelCharacter();
	public char getBoundaryIndelCharacter();
}
