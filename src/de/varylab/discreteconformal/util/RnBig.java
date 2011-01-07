package de.varylab.discreteconformal.util;

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
	
	
	public static BigDecimal[] times(BigDecimal[] dst, BigDecimal[] src1, BigDecimal[] src2) {
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
					BigDecimal prod = src1[i*n+k].multiply(src2[k*n+j]);
					outVal = outVal.add(prod);
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
		BigDecimal[]  src
	)	{
		if (dst == null) dst = new BigDecimal[src.length];
		if (dst.length != src.length) {
			throw new IllegalArgumentException("Vectors must be same length");
		}
		int n = dst.length;
		for (int i=0; i<n; ++i)	{
			dst[i] = src[i].multiply(factor);
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
				if ((x = t[j*n+i].multiply(t[j*n+i])).compareTo(largesq) > 0)	{
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
					t[j*n+k] = t[j*n+k].subtract(f.multiply(t[i*n+k]));  
				}
				for (k=0; k<n; ++k)	{
					minv[j*n+k] = minv[j*n+k].subtract(f.multiply(minv[i*n+k]));
				}
			}
		}

		for (i = 0; i < n; i++) {
			f = t[i*n+i];
			for (f = ONE.divide(f, context), j = 0; j < n; j++) {
				t[i*n+j] = t[i*n+j].multiply(f);
				minv[i*n+j] = minv[i*n+j].multiply(f);
			}
		}
		for (i = n-1; i >= 0; i--)
			for (j = i-1; j >= 0; j--) {
				f = t[j*n+i];
				for (k=0; k<n; ++k)	{
					t[j*n+k] = t[j*n+k].subtract(f.multiply(t[i*n+k]));  
				}
				for (k=0; k<n; ++k)	{
					minv[j*n+k] = minv[j*n+k].subtract(f.multiply(minv[i*n+k]));
				}
			}
		return minv;
	}
	
	
}
