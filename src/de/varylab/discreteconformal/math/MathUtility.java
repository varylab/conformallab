package de.varylab.discreteconformal.math;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.util.logging.Logger;

public class MathUtility {

	private static Logger
		log = Logger.getLogger(MathUtility.class.getName());
	
	public static double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}

	public static double arcosh(double x) {
		if (x < 1.0) {
			log.warning("clamping x value to 1.0, x = " + x);
			x = 1.0;
		}
		double r = x + sqrt(x*x - 1);
		return log(r);
	}
	
}
