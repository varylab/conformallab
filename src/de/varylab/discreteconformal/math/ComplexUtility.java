package de.varylab.discreteconformal.math;

import de.jtem.mfc.field.Complex;

public class ComplexUtility {

	/**
	 * Inverse stereographic projection
	 * @param Z
	 * @param result if given the result s written to this array
	 * @return
	 */
	public static double[] inverseStereographic(Complex Z, double[]... result) {
		if (Z.re == Double.POSITIVE_INFINITY) {
			return new double[] {0,0,1};
		}
		double X = Z.getRe();
		double Y = Z.getIm();
		double d = 1 + X*X + Y*Y;
		double x = 2*X / d;
		double y = 2*Y / d;
		double z = (X*X + Y*Y - 1) / d;
		if (result.length != 0) {
			result[0][0] = x;
			result[0][1] = y;
			result[0][2] = z;
			return result[0];
		}
		return new double[] {x, y, z};
	}

	public static Complex stereographic(double[] p) {
		return new Complex(p[0] / (1 - p[2]), p[1] / (1 - p[2]));
	}

}
