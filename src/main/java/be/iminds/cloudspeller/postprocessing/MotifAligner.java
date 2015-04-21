package be.iminds.cloudspeller.postprocessing;

import java.util.HashMap;
import java.util.Map;

import be.iminds.cloudspeller.alphabets.Alphabet;
import be.iminds.cloudspeller.alphabets.CharacterIterator;


public class MotifAligner {
	
	private Map<String,Double> characterDistanceMap = new HashMap<String,Double>();
	private Alphabet alphabet;
	CharacterDistanceCalculator calculator;
	
	public MotifAligner(Alphabet alphabet, CharacterDistanceCalculator calculator){
		setAlphabet(alphabet);
		this.calculator=calculator;
		preprocessCharacterDistances();
	}

	private void preprocessCharacterDistances() {
		//System.out.println("STart preprocess");
		String allChars =alphabet.getAllChars()+calculator.getBulkIndelCharacter()+calculator.getBoundaryIndelCharacter();
		
		CharacterIterator charIterator1 = new CharacterIterator(allChars);
		while (charIterator1.hasNext()){
			char c1 = charIterator1.next();
			CharacterIterator charIterator2 = new CharacterIterator(allChars);
			while (charIterator2.hasNext()){
				char c2 = charIterator2.next();
				String key = ""+c1+c2;
				//System.out.println(key);
				characterDistanceMap.put(key,calculator.calculateDistance(c1, c2));
			}
		}
		//System.out.println("End preprocess");
	}
	
	public void setAlphabet(Alphabet alphabet){
		this.alphabet=alphabet;
	}
	
	private double minValue3(double v1, double v2, double v3){
		double v12 = Math.min(v1,v2);
		return Math.min(v12,v3);
	}
	
	public double calculateSimilarityDistance(String s1, String s2){
		double [][] grid = generateManhattanGrid(s1,s2);
		processInitialBoundaries(s1,s2,grid);
		processBulk(s1,s2,grid);
		processFinalBoundaries(s1,s2,grid);
		
		int lastRow = grid.length-1;
		int lastCol = grid[lastRow].length-1;
		
		return grid[lastRow][lastCol];
	}

	private void processBulk(String s1, String s2, double [][] grid) {
		
		int secondRow = 1;
		int secondCol = 1;
		int lastRow = grid.length-1;
		int lastCol = grid[lastRow].length-1;
		
		for (int i=secondRow; i<lastRow; i++){
			for (int j=secondCol; j<lastCol; j++){
					
				double diagonalTrack   = grid[i-1][j-1];
				double horizontalTrack = grid[i][j-1];
				double verticalTrack   = grid[i-1][j];
				
				diagonalTrack  +=characterDistanceMap.get(""+s1.charAt(i-1) + s2.charAt(j-1));
				horizontalTrack+=characterDistanceMap.get(""+calculator.getBulkIndelCharacter()+ s2.charAt(j-1));
				verticalTrack  +=characterDistanceMap.get(""+s1.charAt(i-1)+calculator.getBulkIndelCharacter());
				
				grid[i][j] = minValue3(diagonalTrack,horizontalTrack,verticalTrack);
			}
		}
	}

	private void processInitialBoundaries(String s1, String s2, double [][] grid) {
		
		int firstRow = 0;
		int firstCol = 0;
		int lastRow = grid.length-1;
		int lastCol = grid[lastRow].length-1;
		
		//first el
		grid[firstRow][firstCol]=0;
		
		//first column
		for (int i=firstRow+1; i<=lastRow; i++){
			
			grid[i][firstCol] = grid[i-1][firstCol];
			grid[i][firstCol]+=characterDistanceMap.get(""+s1.charAt(i-1)+calculator.getBoundaryIndelCharacter()); 
		}
		
		//first row
		for (int j=firstCol+1; j<=lastCol; j++){
			grid[firstRow][j] = grid[firstRow][j-1];
			grid[firstRow][j]+=characterDistanceMap.get(""+calculator.getBoundaryIndelCharacter()+s2.charAt(j-1));
		}
		
	}
	
	private void processFinalBoundaries(String s1, String s2, double[][] grid) {
		
		int firstRow = 0;
		int firstCol = 0;
		int lastRow = grid.length-1;
		int lastCol = grid[lastRow].length-1;
		
		grid[firstRow][firstCol]=0;
		
		//last column
		for (int i=firstRow+1; i<lastRow; i++){
			
			double diagonalTrack   = grid[i-1][lastCol-1];
			double horizontalTrack = grid[i][lastCol-1];
			double verticalTrack   = grid[i-1][lastCol];
			
			diagonalTrack  +=characterDistanceMap.get(""+s1.charAt(i-1) + s2.charAt(lastCol-1));
			horizontalTrack+=characterDistanceMap.get(""+calculator.getBoundaryIndelCharacter()+ s2.charAt(lastCol-1));
			verticalTrack  +=characterDistanceMap.get(""+s1.charAt(i-1)+calculator.getBoundaryIndelCharacter());
			
			grid[i][lastCol] = minValue3(diagonalTrack,horizontalTrack,verticalTrack);
		}
		
		//last row
		for (int j=firstCol+1; j<lastCol; j++){
			
			double diagonalTrack   = grid[lastRow-1][j-1];
			double horizontalTrack = grid[lastRow][j-1];
			double verticalTrack   = grid[lastRow-1][j];
			
			diagonalTrack  +=characterDistanceMap.get(""+s1.charAt(lastRow-1) + s2.charAt(j-1));
			horizontalTrack+=characterDistanceMap.get(""+calculator.getBoundaryIndelCharacter()+ s2.charAt(j-1));
			verticalTrack  +=characterDistanceMap.get(""+s1.charAt(lastRow-1)+calculator.getBoundaryIndelCharacter());
			
			grid[lastRow][j] = minValue3(diagonalTrack,horizontalTrack,verticalTrack);
		}
		
		//last el
		double diagonalTrack   = grid[lastRow-1][lastCol-1];
		double horizontalTrack = grid[lastRow][lastCol-1];
		double verticalTrack   = grid[lastRow-1][lastCol];
		
		diagonalTrack  +=characterDistanceMap.get(""+s1.charAt(lastRow-1) + s2.charAt(lastCol-1));
		horizontalTrack+=characterDistanceMap.get(""+calculator.getBoundaryIndelCharacter()+ s2.charAt(lastCol-1));
		verticalTrack  +=characterDistanceMap.get(""+s1.charAt(lastRow-1)+calculator.getBoundaryIndelCharacter());
		
		grid[lastRow][lastCol] = minValue3(diagonalTrack,horizontalTrack,verticalTrack);
	}

	private double[][] generateManhattanGrid(String s1, String s2)  {
		int nRows = s1.length()+1;
		int nCols = s2.length()+1;
		double [][] grid = new double [nRows][nCols];
		return grid;
	}
	
}
