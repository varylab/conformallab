package de.varylab.discreteconformal.math;


/**
 * Evaluate Clausens integral
 * <p>
 * Copyright 2005 <a href="http://www.math.tu-berlin.de/~springb/">Boris Springborn</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Boris Springborn
 */
public class Clausen {

	static final double TWO_PI = 2 * Math.PI;
	static final double MINUS_LOG_2 = Math.log(0.5);
	static final double BORDERLINE = 2.0944;
	
	private Clausen() {}

	public static double clausen(double x) {
		x = Math.IEEEremainder(x, TWO_PI);
		if (x == 0.0) {
			return 0.0;
		}
		if (Math.abs(x) <= BORDERLINE) {
			final double xx = x * x;
			return ((((((((((((
					  2.3257441143020875e-22 * xx
					+ 1.0887357368300848e-20) * xx 
					+ 5.178258806090624e-19) * xx
					+ 2.5105444608999545e-17) * xx
					+ 1.2462059912950672e-15) * xx
					+ 6.372636443183181e-14) * xx
					+ 3.387301370953521e-12) * xx
					+ 1.8978869988971e-10) * xx
					+ 1.1482216343327455e-8) * xx
					+ 7.873519778281683e-7) * xx
					+ 0.00006944444444444444) * xx
					+ 0.013888888888888888) * xx 
					- Math.log(Math.abs(x)) + 1.0) * x;
		}
		x += ((x > 0.0) ? - Math.PI : Math.PI);
		final double xx = x * x;
		return ((((((((((((  
				  3.901950904063069e-15 * xx
				+ 4.566487567193635e-14) * xx
				+ 5.429792727596476e-13) * xx
				+ 6.5812165661369675e-12) * xx
				+ 8.167010963952222e-11) * xx
				+ 1.0440290284867003e-9) * xx
				+ 1.3870999114054669e-8) * xx
				+ 1.941538399871733e-7) * xx
				+ 2.927965167548501e-6) * xx
				+ 0.0000496031746031746) * xx
				+ 0.0010416666666666667) * xx
				+ 0.041666666666666664) * xx 
				+ MINUS_LOG_2) * x;
	}
	
}

