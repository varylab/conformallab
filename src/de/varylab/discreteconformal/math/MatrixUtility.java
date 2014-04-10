package de.varylab.discreteconformal.math;

import de.jreality.math.Rn;

public class MatrixUtility {

	public static double[] makeMappingMatrix(double[][] from, double[][] to) {
		double[] from1 = Rn.convertArray2DToArray1D(null, from);
		Rn.transpose(from1, from1);
		Rn.inverse(from1, from1);
		double[] to1 = Rn.convertArray2DToArray1D(null, to);
		Rn.transpose(to1, to1);
		return Rn.times(null, to1, from1);
	}

}
