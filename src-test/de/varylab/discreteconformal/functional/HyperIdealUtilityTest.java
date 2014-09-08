package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.Clausen.Л;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolume;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_13;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static org.junit.Assert.*;

import org.junit.Assert;
import org.junit.Test;

public class HyperIdealUtilityTest {

	@Test
	public void testZeta13() throws Exception {
		double ζ_13 = ζ_13(0.1, 0.1, 0.1);
		double ζ_e = 0.1;
		Assert.assertEquals(ζ_e, ζ_13, 1E-15);
	}
	
	@Test
	public void testVolumeEuclidean() throws Exception {
		double βi = acos(1.0 / 3);
		double βj = acos(1.0 / 3);
		double βk = acos(1.0 / 3);
		double V = calculateTetrahedronVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(0.0, V, 1E-7);
	}
	
	@Test
	public void testVolumeRegularIdeal1() throws Exception {
		double βi = PI / 3;
		double βj = PI / 3;
		double βk = PI / 3;
		double Ve = Л(βi) + Л(βj) + Л(βk);
		double V = calculateTetrahedronVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeRegularIdeal2() throws Exception {
		double βi = PI / 2;
		double βj = PI / 4;
		double βk = PI / 4;
		double Ve = Л(βi) + Л(βj) + Л(βk);
		double V = calculateTetrahedronVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeOctahedron() throws Exception {
		double Ve = 8*Л(PI/4);
		double V = calculateTetrahedronVolume(0, 0, 0, 0, 0, 0);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeHyperidealSingleVertex() throws Exception {
		double βi = PI / 5;
		double βj = PI / 4;
		double βk = PI / 4;
		double ai = PI + βi - βj - βk;
		double aj = PI + βj - βi - βk;
		double ak = PI + βk - βi - βj;
		double aijk = PI - βk - βi - βj;
		double Ve = 0.5 * (Л(βi) + Л(βj) + Л(βk) + Л(ai/2) + Л(aj/2) + Л(ak/2) + Л(aijk/2));
		double V = calculateTetrahedronVolume(βi, βj, βk, ai/2, aj/2, ak/2);
		assertEquals(Ve, V, 1E-12);
	}
	
}
