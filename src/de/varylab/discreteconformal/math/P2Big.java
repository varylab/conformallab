package de.varylab.discreteconformal.math;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

import java.math.BigDecimal;
import java.math.MathContext;

public class P2Big {

	
	public static BigDecimal[] makeDirectIsometryFromFrames(
		BigDecimal[] dst,
		BigDecimal[] p0, 
		BigDecimal[] p1, 
		BigDecimal[] q0, 
		BigDecimal[] q1, 
		int signature,
		MathContext context
	) {
		BigDecimal[] toP = makeDirectIsometryFromFrame(null, p0, p1, signature, context);
		BigDecimal[] toQ = makeDirectIsometryFromFrame(null, q0, q1, signature, context);
		BigDecimal[] iToP = RnBig.inverse(null, toP, context);
		dst = RnBig.times(dst, toQ, iToP, context);
		return dst;
	}

	private static BigDecimal[] makeDirectIsometryFromFrame(
		BigDecimal[] dst, 
		BigDecimal[] p0,
		BigDecimal[] p1, 
		int signature,
		MathContext context
	) {
		if (dst == null) dst = new BigDecimal[9];
		PnBig.normalize(p0, p0, signature, context);
		BigDecimal[] polarP = PnBig.polarize(null, p0, signature);
		BigDecimal[] lineP = lineFromPoints(null, p0, p1, context);
		BigDecimal[] p1n = PnBig.normalize(null, pointFromLines(null, polarP, lineP, context), signature, context);
		BigDecimal[] p2 = PnBig.polarize(null, lineP, signature);
		PnBig.normalize(p2, p2, signature, context);
		makeMatrixFromColumns(dst, p0, p1n, p2);
		return dst;
	}

	
	public static BigDecimal[] pointFromLines(
		BigDecimal[] point,
		BigDecimal[] l1, 
		BigDecimal[] l2,
		MathContext context
	)	{
		if (l1.length < 3 || l2.length < 3)	{
			throw new IllegalArgumentException("Input arrays too short");
		}
		if (point == null) point = new BigDecimal[3];	
		point[0] = l1[1].multiply(l2[2], context).subtract(l1[2].multiply(l2[1], context), context);
		point[1] = l1[2].multiply(l2[0], context).subtract(l1[0].multiply(l2[2], context), context);
		point[2] = l1[0].multiply(l2[1], context).subtract(l1[1].multiply(l2[0], context), context);
		return point;
	}
	
	/**
	 * Calculate the line coordinates of the line connecting the two points <i>p1</i> and <i>p2</i>.
	 * @param point
	 * @param l1
	 * @param l2
	 * @return
	 */
	public static BigDecimal[] lineFromPoints(BigDecimal[] line, BigDecimal[] p1, BigDecimal[] p2, MathContext context)	{
		return pointFromLines(line, p1, p2, context);
	}
	
	
	public static BigDecimal[] makeMatrixFromColumns(BigDecimal[] dst, BigDecimal[] p0, BigDecimal[] p1, BigDecimal[] p2) {
		if (dst == null) dst = new BigDecimal[9];
		BigDecimal[][] ptrs = {p0, p1, p2};
		for (int i = 0; i<3; ++i)	{
			for (int j = 0; j<3; ++j)	{
				dst[3*i+j] = ptrs[j][i];
			}
		}
		return dst;
	}
	
	
	public static BigDecimal[] projectP3ToP2(BigDecimal[] vec3, BigDecimal[] vec4)	{
		BigDecimal[] dst;
		if (vec3 == null) {
			dst = new BigDecimal[3];
		}
		else dst = vec3;
		dst[0] = vec4[0];
		dst[1] = vec4[1];
		dst[2] = vec4[3];
		return dst;
	}
	
	/**
	 * Convert (x,y,z) into (x,y,0,z)
	 * @param vec4
	 * @param vec3
	 * @return
	 */
	 public static BigDecimal[] imbedP2InP3(BigDecimal[] vec4, BigDecimal[] vec3)	{
		BigDecimal[] dst;
		if (vec4 == null)	{
			dst = new BigDecimal[4];
		}
		else dst = vec4;
		dst[0] = vec3[0];
		dst[1] = vec3[1];
		dst[3] = vec3[2];
		dst[2] = ZERO;
		return dst;
	}

	 private static int[] which = {0,1,3};
	 public static BigDecimal[] imbedMatrixP2InP3(BigDecimal[] dst, BigDecimal[] m3)	{
		if (dst == null)	dst = new BigDecimal[16];
		for (int i = 0; i<3; ++i)	{
			int i4 = which[i];
			for (int j = 0; j<3; ++j)	{
				int j4 = which[j] + 4 * i4;
				int j3 = i*3 + j;
				dst[j4] = m3[j3];
			}
		}
		dst[2] = dst[6] = dst[8] = dst[9] = dst[11] = dst[14] = ZERO; dst[10] = ONE;
		return dst;
	}
	
	
}
