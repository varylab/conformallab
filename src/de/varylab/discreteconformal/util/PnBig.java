package de.varylab.discreteconformal.util;

import static de.varylab.discreteconformal.util.RnBig.isOne;
import static de.varylab.discreteconformal.util.RnBig.isZero;
import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

import java.math.BigDecimal;
import java.math.MathContext;

import de.jreality.math.Pn;

public class PnBig {

	public static final int ELLIPTIC = 1;
	public static final int EUCLIDEAN = 0;
	public static final int HYPERBOLIC = -1;
	public static final int PROJECTIVE = 2;
	public static BigDecimal[] zDirectionP3 = { 
		new BigDecimal(0.0),
		new BigDecimal(0.0),
		new BigDecimal(1.0), 
		new BigDecimal(0.0) 
	};
	
	public static BigDecimal[] dehomogenize(
		BigDecimal[] dst, 
		BigDecimal[] src,
		MathContext context
	)	{
		// assert dim checks
		int sl = src.length;
		if (dst == null) dst = new BigDecimal[src.length];
		int dl = dst.length;
		// allow dst array same length or one shorter than source
		if (! (dl == sl || dl +1 == sl))	{
			throw new IllegalArgumentException("Invalid dimensions");
		}
		BigDecimal last = src[sl-1];
		if (isOne(last) || isZero(last)) 	{
			if (src != dst) {
				System.arraycopy(src,0,dst,0,dl);
			}
			return dst;
		}
		last = ONE.divide(last, context);
		for (int i = 0; i<dl; ++i)		{
			dst[i] = last.multiply(src[i]);
		}
		if (dl == sl) {
			dst[dl-1] = ONE;
		}
		return dst;
	}
	
	public static BigDecimal norm(BigDecimal[] src, int metric, MathContext context) {
		BigDecimal d = innerProduct(src, src, metric, context).abs();
		double sqrt = Math.sqrt(d.doubleValue());
		return new BigDecimal(sqrt);
	}
	
	public static BigDecimal innerProduct(
		BigDecimal dst[], 
		BigDecimal[] src, 
		int metric,
		MathContext context
	)	{
		// assert dim checks
		BigDecimal sum = new BigDecimal(0.0);
		if (src.length != dst.length)		{
			throw new IllegalArgumentException("Incompatible lengths");
		}
		int n = dst.length;
		
		for (int i = 0; i< n-1; ++i) {
			sum = sum.add(dst[i].multiply(src[i], context), context);
		}
		BigDecimal ff = dst[n-1].multiply(src[n-1], context);
		
		switch(metric)	{
			case HYPERBOLIC:		// metric (n-1,1)
				sum = sum.subtract(ff, context);
				break;
			case EUCLIDEAN:		// metric (n-1, 0)				
				if (!(RnBig.isOne(ff) || RnBig.isZero(ff)))	{
					sum = sum.divide(ff, context);
				}
				break;
			case ELLIPTIC:		// metric (n, 0)
				sum = sum.add(ff, context);
				break;
		}
		return sum;
	}
	

	public static BigDecimal[] normalize(
		BigDecimal[] dst, 
		BigDecimal[] src,
		int metric,
		MathContext context
	) {
		if (dst == null)
			dst = new BigDecimal[src.length];
		if (metric == EUCLIDEAN) {
			return dehomogenize(dst, src, context);
		}
		return setToLength(dst, src, ONE, metric, context);
	}

	public static BigDecimal[] setToLength(
		BigDecimal[] dst, 
		BigDecimal[] src,
		BigDecimal length, 
		int metric,
		MathContext context
	) {
		// assert dim checks
		if (dst == null)
			dst = new BigDecimal[src.length];
		if (dst.length != src.length) {
			throw new IllegalArgumentException("Incompatible lengths");
		}
		if (metric == EUCLIDEAN) {
			dehomogenize(dst, src, context);
		} else {
			System.arraycopy(src, 0, dst, 0, dst.length);
		}
		BigDecimal ll = norm(dst, metric, context);
		if (ll.compareTo(ZERO) == 0) {
			return dst;
		}
		ll = length.divide(ll, context);
		RnBig.times(dst, ll, dst, context);
		if (metric == EUCLIDEAN && dst[dst.length - 1].compareTo(ZERO) != 0) {
			dst[dst.length - 1] = new BigDecimal(1.0);
		}
		return dst;
	}

	public static BigDecimal[] polarize(BigDecimal[] polar, BigDecimal[] p, int metric)	{
		if (polar == null)	polar = (BigDecimal[]) p.clone();
		else System.arraycopy(p,0,polar, 0, p.length);
		// last element is multiplied by the metric!
		switch (metric)	{
			case Pn.ELLIPTIC:
				// self-polar!
				break;
			case Pn.EUCLIDEAN:
				polar[polar.length - 1] = ZERO;
				break;
			case Pn.HYPERBOLIC:
				polar[polar.length-1] = polar[polar.length-1].negate();
				break;
		}
		return polar;
	}
	
	
}
