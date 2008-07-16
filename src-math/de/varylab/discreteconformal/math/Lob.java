package de.varylab.discreteconformal.math;

public class Lob {

	public static double valueAt(double x) {
		return 2 * Clausen.valueAt(x / 2);
	}
	
}
