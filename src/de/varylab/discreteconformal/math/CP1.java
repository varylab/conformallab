package de.varylab.discreteconformal.math;

import de.jreality.math.Rn;
import de.jtem.mfc.field.Complex;

public class CP1 {

	protected static double[] sphereToParaboloid = {
		1, 0, 0, 1,
		0, 1, 0, 0, 
		0, 0, 1, 0, 
		-1, 0, 0, 1
	};
	
	protected static double[] paraboloidToSphere = {
		1, 0, 0, -1,
		0, 2, 0, 0, 
		0, 0, 2, 0, 
		1, 0, 0, 1
	};
	
	protected static double[] zxyToxyz = {
		0,1,0,0,
		0,0,1,0,
		1,0,0,0,
		0,0,0,1
	};
	
	public static double[] convertPSL2CToSO31(double[] dst, Complex[] lft)	{
		if (dst == null) {
			dst = new double[16];
		}
		double ax = lft[0].re;
		double ay = lft[0].im;
		double bx = lft[1].re;
		double by = lft[1].im;
		double cx = lft[2].re;
		double cy = lft[2].im;
		double dx = lft[3].re;
		double dy = lft[3].im;
		double[] tmp = new double[16];
		tmp[0] = ax*ax + ay*ay;
		tmp[1] = 2*(ax*bx+ay*by);
		tmp[2] = 2*(ax*by-ay*bx);
		tmp[3] = bx*bx+by*by;
		tmp[4] = ax*cx+ay*cy;
		tmp[5] = bx*cx + by*cy + ax*dx + ay*dy;
		tmp[6] = by*cx - bx*cy - ay*dx + ax*dy;
		tmp[7] = bx*dx + by*dy;
		tmp[8] = ay*cx - ax*cy;
		tmp[9] = by*cx - bx*cy + ay*dx - ax*dy;
		tmp[10]= -(bx*cx) - by*cy + ax*dx + ay*dy;
		tmp[11]= by*dx - bx*dy;
		tmp[12]= cx*cx+cy*cy;
		tmp[13]= 2*(cx*dx + cy*dy);
		tmp[14]= 2*(-cy*dx + cx*dy);
		tmp[15]= dx*dx+dy*dy;
		Rn.conjugateByMatrix(dst, tmp, Rn.times(null, zxyToxyz, paraboloidToSphere));
		return dst;
	}
	
	public static Complex[] projectivity(Complex[] dst, Complex z1, Complex z2, Complex z3, Complex w1, Complex w2, Complex w3)	{
		//TODO make sure inputs are valid
		Complex[] m1 = standardProjectivity(null, z1, z2, z3);
		Complex[] m2 = standardProjectivity(null, w1, w2, w3);
		Complex[] im2 = Cn.invert(null, m2);
		dst = Cn.times(dst, im2, m1);
		return dst;
	}
	
	/**
	 * Generate the Moebius tform taking (z1,z2,z3) to (0, 1, infinity)
	 * @param m
	 * @param z1
	 * @param z2
	 * @param z3
	 * @return
	 */
	public static Complex[] standardProjectivity(Complex[] m, Complex z1, Complex z2, Complex z3)	{
		if (m == null || m.length != 4) {
			m = new Complex[4];
		}
		for (int i = 0; i<4; ++i) {
			if (m[i] == null) m[i] = new Complex();
		}
		if (z1.isInfinite())	{
			m[0].assign(0,0);
			m[1].assignMinus(z2, z3);
			m[2].assign(1,0);
			m[3].assignMinus(z3);
		} else if (z2.isInfinite())	{
			m[0].assign(1,0);
			m[1].assignMinus(z1);
			m[2].assign(1,0);
			m[3].assignMinus(z3);
		} else if (z3.isInfinite())	{
			m[0].assign(1,0);
			m[1].assignMinus(z1);
			m[2].assign(0,0);
			m[3].assignMinus(z2, z1);
		} else {
			m[0].assignMinus(z2, z3);
			m[1].assignTimes(z1.neg(), m[0]);
			m[2].assignMinus(z2, z1);
			m[3].assignTimes(z3.neg(), m[2]);
		}
//		Cn.normalize(m,m);
		return m;
	}
	
}
