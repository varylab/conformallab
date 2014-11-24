package de.varylab.discreteconformal.util;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.mfc.field.Complex;

public class DiscreteEllipticUtilityTest {

	@Test
	public void testNormalizeModulus1() throws Exception {
		Complex tau = new Complex(0.45, 1.1);
		Complex tauNormalized = DiscreteEllipticUtility.normalizeModulus(tau);
		Assert.assertEquals(0.3185840707964602, tauNormalized.re, 1E-12);
		Assert.assertEquals(0.7787610619469026, tauNormalized.im, 1E-12);
	}
	
	
	@Test
	public void testNormalizeModulusPeriodShift() throws Exception {
		Complex tau1 = new Complex(0.3, 1.0);
		Complex tau2 = new Complex(-0.7, 1.0);
		Complex tau1Normalized = DiscreteEllipticUtility.normalizeModulus(tau1);
		Complex tau2Normalized = DiscreteEllipticUtility.normalizeModulus(tau2);
		Assert.assertEquals("the moduli should be equal", tau1Normalized, tau2Normalized);
	}
	
}
