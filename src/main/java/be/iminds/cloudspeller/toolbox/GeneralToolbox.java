package be.iminds.cloudspeller.toolbox;

import java.util.Locale;
import java.util.Scanner;

public class GeneralToolbox {

	public static Scanner generateScanner(String s){
		Scanner scanner = new Scanner(s);
		scanner.useLocale(Locale.US); //interprete . as decimal sign -> standard is , !!
		return scanner;
	}
	
	public static void parseConfidenceGraphValues(String line, int [] F, double[] p){

		if (line == null){
			System.err.println("Empy line");
			return;
		}

		Scanner scanner = generateScanner(line);
		
		scanner.next(); scanner.next();
		
		for (int i=0; i<F.length; i++){
			F[i]=Integer.parseInt(scanner.next());
		}
		
		for (int i=0; i<p.length; i++){
			p[i]=Double.parseDouble(scanner.next());
		}
		
		
		
	}
	
}
