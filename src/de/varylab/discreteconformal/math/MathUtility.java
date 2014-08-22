package de.varylab.discreteconformal.math;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

public class MathUtility {

	public static double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}

	public static double arcosh(double x) {
		double r = x + sqrt(x*x - 1);
		return log(r);
	}
	
}
