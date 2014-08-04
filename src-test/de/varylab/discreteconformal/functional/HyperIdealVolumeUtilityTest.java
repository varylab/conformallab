package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.Clausen.Л;
import static de.varylab.discreteconformal.functional.HyperIdealVolumeUtility.calculateVolume;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class HyperIdealVolumeUtilityTest {

	@Test
	public void testEuclidean() throws Exception {
		double βi = acos(1.0 / 3);
		double βj = acos(1.0 / 3);
		double βk = acos(1.0 / 3);
		double V = calculateVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(0.0, V, 1E-7);
	}
	
	@Test
	public void testRegularIdeal1() throws Exception {
		double βi = PI / 3;
		double βj = PI / 3;
		double βk = PI / 3;
		double Ve = Л(βi) + Л(βj) + Л(βk);
		double V = calculateVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testRegularIdeal2() throws Exception {
		double βi = PI / 2;
		double βj = PI / 4;
		double βk = PI / 4;
		double Ve = Л(βi) + Л(βj) + Л(βk);
		double V = calculateVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testOctahedron() throws Exception {
		double Ve = 8*Л(PI/4);
		double V = calculateVolume(0, 0, 0, 0, 0, 0);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testHyperidealSingleVertex() throws Exception {
		double βi = PI / 2;
		double βj = PI / 4;
		double βk = PI / 4;
		double ai = PI + βi - βj - βk;
		double aj = PI + βj - βi - βk;
		double ak = PI + βk - βi - βj;
		double aijk = PI - βk - βi - βj;
		double Ve = 0.5 * (Л(βi) + Л(βj) + Л(βk) + Л(ai/2) + Л(aj/2) + Л(ak/2) + Л(aijk/2));
		double V = calculateVolume(βi, βj, βk, 0, 0, 0);
		assertEquals(Ve, V, 1E-12);
	}
	
}
