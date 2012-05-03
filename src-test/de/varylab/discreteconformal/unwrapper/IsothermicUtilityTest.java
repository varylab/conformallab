package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;

import org.junit.Assert;
import org.junit.Test;

import de.varylab.discreteconformal.unwrapper.IsothermicUtility;

public class IsothermicUtilityTest {

	private final double
		EPS = 1E-2;
	
	@Test
	public void testCalculateTriangleAngle() throws Exception {
		double a1 = IsothermicUtility.calculateTriangleAngle(PI/2, 0, PI/4);
		Assert.assertEquals(PI/2, a1, 1E-10);
		double a2 = IsothermicUtility.calculateTriangleAngle(PI/2, 0, -PI/4);
		Assert.assertEquals(PI/2, a2, 1E-10);
		double a3 = IsothermicUtility.calculateTriangleAngle(-PI/2 + EPS, PI/2 - EPS, 0);
		Assert.assertEquals(2*EPS, a3, 1E-10);
		double a4 = IsothermicUtility.calculateTriangleAngle(-PI/2 + EPS, PI/2 - EPS, PI/2);
		Assert.assertEquals(PI - 2*EPS, a4, 1E-10);
		double a5 = IsothermicUtility.calculateTriangleAngle(-PI/8, -3*PI/8, PI/2);
		Assert.assertEquals(PI/4, a5, 1E-10);
		double a6 = IsothermicUtility.calculateTriangleAngle(PI/8, 3*PI/8, PI/2);
		Assert.assertEquals(PI/4, a6, 1E-10);
		double a7 = IsothermicUtility.calculateTriangleAngle(-PI/4, PI/4, PI/2);
		Assert.assertEquals(PI/2, a7, 1E-10);
		double a8 = IsothermicUtility.calculateTriangleAngle(-PI/4, PI/4, 0);
		Assert.assertEquals(PI/2, a8, 1E-10);
	}
	
}
