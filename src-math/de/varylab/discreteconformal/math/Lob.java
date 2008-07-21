package de.varylab.discreteconformal.math;

public class Lob {

	public static double lob(double x) {
		return 2 * Clausen.clausen(x / 2);
	}
	
}
