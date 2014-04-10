package de.varylab.discreteconformal.math;

import de.jtem.mfc.field.Complex;

public class Cn {

	public static Complex determinant(Complex[] m)	{
		return m[0].times(m[3]).minus(m[1].times(m[2]));
	}
	
	public static Complex[] invert(Complex[] dst, Complex[] src)	{
		if (dst == null) {
			dst = new Complex[4];
			for (int i =0; i<4; ++i) {
				dst[i] = new Complex();
			}
		}
		Complex[] target = dst;
		if (src == dst)		{
			target = new Complex[4];
			for (int i =0; i<4; ++i) {
				target[i] = new Complex();
			}
		}
		target[0].assign(src[3]);
		target[1].assignNeg(src[1]);
		target[2].assignNeg(src[2]);
		target[3].assign(src[0]);
		if (target != dst)	{
			for (int i = 0; i<4; ++i)  {
				dst[i].assign(target[i]);
			}
		}
		return dst;
	}
	
	public static Complex[] times(Complex[] dst, final Complex[] src1, final Complex[] src2)		{
		if (dst == null) {
			dst = new Complex[4];
			for (int i =0; i<4; ++i) dst[i] = new Complex();
		}
		Complex[] target = dst;
		if (src1 == dst || src2 == dst)		{
			target = new Complex[4];
			for (int i =0; i<4; ++i) target[i] = new Complex();
		}
		for (int i = 0; i<2; ++i)	{
			for (int j = 0; j<2; ++j)	{
				target[2*i+j].assign(0, 0);
				for (int k = 0; k<2; ++k) {
					target[i*2+j].assignPlus(target[i*2+j], src1[i*2+k].times(src2[2*k+j]));
				}
			}
		}
		if (target != dst) {
			for (int i = 0; i<4; ++i) {
				dst[i].assign(target[i]);
			}
		}
		return dst;
	}
	
	public static Complex[] times(Complex[] dst, Complex d, Complex[] src) {
		if (dst == null) {
			dst = new Complex[4];
			for (int i =0; i<4; ++i) dst[i] = new Complex();
		}
		int n= src.length;
		for (int i = 0; i<n; ++i)	{
			dst[i].assignTimes(d, src[i]);
		}
		return dst;
	}	

	public static Complex[] normalize(Complex[] dst, Complex[] src)	{
		Complex d = Cn.determinant(src);
		if (d.absSqr() == 0.0) {
			return dst;
		}
		d = d.invert();
		d = d.sqrt();
		dst = times(dst, d, src);
		return dst;
	}

}
