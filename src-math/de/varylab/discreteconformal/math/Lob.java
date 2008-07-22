package de.varylab.discreteconformal.math;

import static java.lang.Math.PI;

public class Lob {

	public static double lob(double x) {
		return (0 <= x && x <= PI) ? Clausen.clausen2(2*x)/2 : -1e10;
	}
	
}
