package de.varylab.discreteconformal.math;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

import java.math.BigDecimal;
import java.math.MathContext;

import de.jreality.math.Matrix;
import de.jreality.math.Rn;

public class RnBig {
	
	public static boolean isZero(BigDecimal d) {
		return d.compareTo(BigDecimal.ZERO) == 0;
	}
	
	public static boolean isOne(BigDecimal d) {
		return d.compareTo(BigDecimal.ONE) == 0;
	}
	
	public static BigDecimal[] toBig(BigDecimal[] dst, Matrix T) {
		return toBig(dst, T.getArray());
	}
	
	public static BigDecimal[] toBig(BigDecimal[] dst, double[] T) {
		if (dst == null || dst.length != T.length) {
			dst = new BigDecimal[T.length];
		}
		for (int i = 0; i < T.length; i++) {
			dst[i] = new BigDecimal(T[i]);
		}
		return dst;
	}
	
	
	public static double[] toDouble(double[] dst, BigDecimal[] T) {
		if (dst == null || dst.length != T.length) {
			dst = new double[T.length];
		}
		for (int i = 0; i < T.length; i++) {
			dst[i] = T[i].doubleValue();
		}
		return dst;
	}
	
	
	public static double[][] toDouble(double[][] dst, BigDecimal[][] T) {
		if (dst == null || dst.length != T.length) {
			dst = new double[T.length][];
		}
		for (int i = 0; i < T.length; i++) {
			if (dst[i] == null || dst[i].length != T[i].length) {
				dst[i] = new double[T[i].length];
			}
			for (int j = 0; j < T[i].length; j++) {
				dst[i][j] = T[i][j].doubleValue();
			}
		}
		return dst;
	}
	
	
	
	public static String matrixToString(BigDecimal[] m) {
		StringBuffer sb = new StringBuffer();
		int n = Rn.mysqrt(m.length);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
//	        sb.append(df.format(v[4*i+j])).append(j == 3 ? "\n":"\t");
				sb.append(String.format("%g", new Object[]{m[n*i+j]}));
				sb.append(j == (n-1) ? "\n":"\t");
			}
		} 
		return sb.toString();
	}
	
	public static BigDecimal[] matrixTimesVector(
		BigDecimal[] dst,
		BigDecimal[] m, 
		BigDecimal[] src, 
		MathContext context
	) {
		// assert dim check
		BigDecimal[] out;
		boolean rewrite = false;
		if (dst == m || dst == src) {
			out = new BigDecimal[dst.length];
			rewrite = true;
		} else if (dst == null) {
			out = new BigDecimal[src.length];
		} else {
			out = dst;
		}
		_matrixTimesVectorSafe(out, m, src, context);
		if (rewrite) {
			System.arraycopy(out, 0, dst, 0, out.length);
			return dst;
		}
		return out;
	}

	/**
	 * This behaves correctly even if <i>src=dst</i>.
	 * 
	 * @param dst
	 * @param m
	 * @param src
	 */
	private static void _matrixTimesVectorSafe(
		BigDecimal[] dst,
		BigDecimal[] m, 
		BigDecimal[] src, 
		MathContext context
	) {
		int sl = src.length;
		int ml = Rn.mysqrt(m.length);
		boolean dehomog = false;
		if (ml == sl + 1)
			dehomog = true;
		if (sl + 1 < ml || sl > ml) {
			throw new IllegalArgumentException(
					"Invalid dimension in _matrixTimesVectorSafe");
		}
		BigDecimal[] out;
		if (dehomog) {
			out = new BigDecimal[ml];
		} else
			out = dst;
		for (int i = 0; i < ml; ++i) {
			out[i] = ZERO;
			for (int j = 0; j < ml; ++j) {
				if (dehomog && j == ml - 1) {
					out[i] = out[i].add(m[i * ml + j], context);
				} else {
					out[i] = out[i].add(m[i * ml + j].multiply(src[j], context), context);
				}
			}
		}
		if (dehomog) {
			PnBig.dehomogenize(dst, out, context);
		}
	}
	
	
	
	public static BigDecimal[] times(BigDecimal[] dst, BigDecimal[] src1, BigDecimal[] src2, MathContext context) {
		if (src1.length != src2.length) {
			throw new IllegalArgumentException("Matrices must be same size");
		}
		int n = Rn.mysqrt(src1.length);
		BigDecimal[] out;
		boolean rewrite = false;
		if (dst == src1 || dst == src2 || dst == null )	{
			out = new BigDecimal[src1.length];
			if (dst != null) rewrite = true;
		}
		else		{
			out = dst;
		}
		if (out.length != src1.length) {
			throw new IllegalArgumentException("Matrices must be same size");
		}
	
		for (int i=0; i<n; ++i)	{	
			for (int j=0; j<n; ++j)	{
				BigDecimal outVal = new BigDecimal(0.0);
				for (int k=0; k<n; ++k)		{
					// the (i,j)th position is the inner product of the ith row and 
					// the jth column of the two factors
					BigDecimal prod = src1[i*n+k].multiply(src2[k*n+j], context);
					outVal = outVal.add(prod, context);
				}
				out[i*n+j] = outVal;
			}
		}
		if (dst == null)	return out;
		if (rewrite) System.arraycopy(out, 0, dst, 0, dst.length);
		return out;
	}
	
	
	public static BigDecimal[] times(
		BigDecimal[]  dst, 
		BigDecimal factor, 
		BigDecimal[]  src,
		MathContext context
	)	{
		if (dst == null) dst = new BigDecimal[src.length];
		if (dst.length != src.length) {
			throw new IllegalArgumentException("Vectors must be same length");
		}
		int n = dst.length;
		for (int i=0; i<n; ++i)	{
			dst[i] = src[i].multiply(factor, context);
		}
		return dst;
	}
	
	
	public static BigDecimal[] setIdentityMatrix(BigDecimal[] mat) {
		int n = Rn.mysqrt(mat.length), noffs = n + 1;
		for (int i = 0; i < mat.length; i++) {
			mat[i] = ZERO;
		}
		for (int i = 0, k = 0; i < n; i++, k += noffs) {
			mat[k] = ONE;
		}
		return mat;
	}
	
	public static BigDecimal[] inverse(
		BigDecimal[] minvIn, 
		BigDecimal[] m,
		MathContext context
	)	{
		int n = Rn.mysqrt(m.length);
		int i, j, k;
		BigDecimal x = ZERO, f = ZERO;
		BigDecimal[] t;
		int largest;

		t =  new BigDecimal[m.length];
		System.arraycopy(m,0,t,0,m.length);
		
		BigDecimal[] minv;
		if (minvIn == null)	{
			minv = new BigDecimal[m.length];
		} else {
			minv = minvIn;
		}
		setIdentityMatrix(minv);
		
		for (i = 0; i < n; i++) {
			largest = i;
			BigDecimal largesq = t[n*i+i].multiply(t[n*i+i]);
			// find the largest entry in the ith column below the ith row
			for (j = i+1; j < n; j++)
				if ((x = t[j*n+i].multiply(t[j*n+i], context)).compareTo(largesq) > 0)	{
					largest = j;  
					largesq = x;
				}

			for (k=0; k<n; ++k)		{		// swap the ith row with the row with largest entry in the ith column
				x = t[i*n+k]; 
				t[i*n+k] = t[largest*n+k]; 
				t[largest*n+k] = x;
			}
			for (k=0; k<n; ++k)		{		// do the same in the inverse matrix
				x = minv[i*n+k]; 
				minv[i*n+k] = minv[largest*n+k]; 
				minv[largest*n+k] = x;
			}
			// now for each remaining row, subtract off a multiple of the ith row so that the
			// entry in the ith column of that row becomes 0.  Previous entries are already 0
			for (j = i+1; j < n; j++) {
				f = t[j*n+i].divide(t[i*n+i], context);
				for (k=0; k<n; ++k)	{
					t[j*n+k] = t[j*n+k].subtract(f.multiply(t[i*n+k]), context);  
				}
				for (k=0; k<n; ++k)	{
					minv[j*n+k] = minv[j*n+k].subtract(f.multiply(minv[i*n+k]), context);
				}
			}
		}

		for (i = 0; i < n; i++) {
			f = t[i*n+i];
			for (f = ONE.divide(f, context), j = 0; j < n; j++) {
				t[i*n+j] = t[i*n+j].multiply(f, context);
				minv[i*n+j] = minv[i*n+j].multiply(f, context);
			}
		}
		for (i = n-1; i >= 0; i--)
			for (j = i-1; j >= 0; j--) {
				f = t[j*n+i];
				for (k=0; k<n; ++k)	{
					t[j*n+k] = t[j*n+k].subtract(f.multiply(t[i*n+k]), context);  
				}
				for (k=0; k<n; ++k)	{
					minv[j*n+k] = minv[j*n+k].subtract(f.multiply(minv[i*n+k], context), context);
				}
			}
		return minv;
	}
	
	
	public static BigDecimal[] linearCombination(
		BigDecimal[] dst,
		BigDecimal a,
		BigDecimal[] aVec,
		BigDecimal b, 
		BigDecimal[] bVec,
		MathContext context
	) {
		if (aVec.length != bVec.length) {
			throw new IllegalArgumentException("Vectors must be same length");
		}
		if (dst == null) {
			dst = new BigDecimal[aVec.length];
		}
		BigDecimal[] tmp = new BigDecimal[dst.length];
		return add(dst, times(tmp, a, aVec, context),
				times(dst, b, bVec, context), context);
	}

	public static BigDecimal[] add(
		BigDecimal[] dst, 
		BigDecimal[] src1,
		BigDecimal[] src2, 
		MathContext context
	) {
		if (dst == null) {
			dst = new BigDecimal[src1.length];
		}
		int n = src1.length;
		if (src1.length != src2.length) {
			n = Math.min(Math.min(dst.length, src1.length), src2.length);
		}
		for (int i = 0; i < n; ++i) {
			dst[i] = src1[i].add(src2[i], context);
		}
		return dst;
	}
	
}
