package de.varylab.discreteconformal.math;

import java.math.BigDecimal;
import java.math.MathContext;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.P2;
import de.jreality.math.Pn;

public class P2BigTest {

	public static MathContext 
		context = new MathContext(50);
	
	public BigDecimal b(double v) {
		return new BigDecimal(v);
	}
	
	@Test
	public void testMakeDirectIsometryFromFramesEuclidean() throws Exception {
		double[] s1 = {-1.4142135623730963, 0.0, 1.0};
		double[] s2 = {1.4142135623730951, 0.0, 1.0};
		double[] t1 = {-2.828427124746189, 2.4494897427831805, 1.0};
		double[] t2 = {0.0, 2.4494897427831783, 1.0};
		double[] T = P2.makeDirectIsometryFromFrames(null, s1, s2, t1, t2, Pn.EUCLIDEAN);
		
		BigDecimal[] s1Big = {b(-1.4142135623730963), b(0.0), b(1.0)};
		BigDecimal[] s2Big = {b(1.4142135623730951), b(0.0), b(1.0)};
		BigDecimal[] t1Big = {b(-2.828427124746189), b(2.4494897427831805), b(1.0)};
		BigDecimal[] t2Big = {b(0.0), b(2.4494897427831783), b(1.0)};
		BigDecimal[] TBig = P2Big.makeDirectIsometryFromFrames(null, s1Big, s2Big, t1Big, t2Big, Pn.EUCLIDEAN, context);
		
		for (int i = 0; i < TBig.length; i++) {
			Assert.assertEquals(T[i], TBig[i].doubleValue(), 1E-10);
		}
	}
	
	
	
}
