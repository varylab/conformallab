package de.varylab.discreteconformal.plugin;

import static java.lang.Math.cosh;
import static java.lang.Math.sinh;

import org.junit.Assert;
import org.junit.Test;

public class HyperIdealVisualizationPluginTest extends Assert {

	@Test
	public void testGetEuclideanCircleFromHyperbolic_Centered() throws Exception {
		double[] center = {0.0, 0.0, 0.0, 1.0};
		double radius = 1.0;
		double[] result = HyperIdealVisualizationPlugin.getEuclideanCircleFromHyperbolic(center, radius);
		assertEquals(0.0, result[0], 1E-12);
		assertEquals(0.0, result[1], 1E-12);
		assertEquals(sinh(1.0) / (cosh(1.0) + 1), result[2], 1E-12);
	}
	
	@Test
	public void testGetEuclideanCircleFromHyperbolic_OffCenter() throws Exception {
		double[] center = {sinh(1.0), 0.0, 0.0, cosh(1.0)};
		double radius = 1.0;
		double[] result = HyperIdealVisualizationPlugin.getEuclideanCircleFromHyperbolic(center, radius);
		assertEquals(result[2], result[0], 1E-12);
		assertEquals(0.0, result[1], 1E-12);
	}
	
}
